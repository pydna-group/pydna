# -*- coding: utf-8 -*-
"""
This module provides classes that roughly map to the `OpenCloning <https://opencloning.org>`_
data model, which is defined using `LinkML <https://linkml.io>`, and available as a python
package `opencloning-linkml <https://pypi.org/project/opencloning-linkml/>`_. These classes
are documented there, and the ones in this module essentially replace the fields pointing to
sequences and primers (which use ids in the data model) to ``Dseqrecord`` and ``Primer``
objects, respectively. Similarly, it uses Location from ``Biopython`` instead of a string,
which is what the data model uses.

When using pydna to plan cloning, it stores the provenance of ``Dseqrecord`` objects in
their ``source`` attribute. Not all methods generate sources so far, so refer to the
documentation notebooks for examples on how to use this feature. The ``history`` method of
``Dseqrecord`` objects can be used to get a string representation of the provenance of the
sequence. You can also use the ``CloningStrategy`` class to create a JSON representation of
the cloning strategy. That ``CloningStrategy`` can be loaded in the OpenCloning web interface
to see a representation of the cloning strategy.


Contributing
============

Not all fields can be readily serialized to be converted to regular types in pydantic. For
instance, the ``coordinates`` field of the ``GenomeCoordinatesSource`` class is a
``SimpleLocation`` object, or the ``input`` field of ``Source`` is a list of ``SourceInput``
objects, which can be ``Dseqrecord`` or ``Primer`` objects, or ``AssemblyFragment`` objects.
For these type of fields, you have to define a ``field_serializer`` method to serialize them
to the correct type.

"""

from __future__ import annotations

from typing import Optional, Union, Any, ClassVar, Type
from pydantic_core import core_schema
from contextlib import contextmanager
from threading import local

from pydantic import BaseModel, ConfigDict, Field, field_serializer, field_validator

from opencloning_linkml.datamodel import (
    CloningStrategy as _BaseCloningStrategy,
    DatabaseSource as _DatabaseSource,
    Primer as _PrimerModel,
    Source as _Source,
    TextFileSequence as _TextFileSequence,
    AssemblySource as _AssemblySource,
    SourceInput as _SourceInput,
    AssemblyFragment as _AssemblyFragment,
    ManuallyTypedSource as _ManuallyTypedSource,
    RestrictionAndLigationSource as _RestrictionAndLigationSource,
    GibsonAssemblySource as _GibsonAssemblySource,
    RestrictionEnzymeDigestionSource as _RestrictionEnzymeDigestionSource,
    SequenceCutSource as _SequenceCutSource,
    RestrictionSequenceCut as _RestrictionSequenceCut,
    SequenceCut as _SequenceCut,
    InFusionSource as _InFusionSource,
    OverlapExtensionPCRLigationSource as _OverlapExtensionPCRLigationSource,
    InVivoAssemblySource as _InVivoAssemblySource,
    LigationSource as _LigationSource,
    GatewaySource as _GatewaySource,
    GatewayReactionType,
    AnnotationTool,
    HomologousRecombinationSource as _HomologousRecombinationSource,
    CreLoxRecombinationSource as _CreLoxRecombinationSource,
    PCRSource as _PCRSource,
    CRISPRSource as _CRISPRSource,
    RepositoryIdSource as _RepositoryIdSource,
    UploadedFileSource as _UploadedFileSource,
    AddgeneIdSource as _AddgeneIdSource,
    AddgeneSequenceType,
    BenchlingUrlSource as _BenchlingUrlSource,
    SnapGenePlasmidSource as _SnapGenePlasmidSource,
    EuroscarfSource as _EuroscarfSource,
    WekWikGeneIdSource as _WekWikGeneIdSource,
    SEVASource as _SEVASource,
    IGEMSource as _IGEMSource,
    OpenDNACollectionsSource as _OpenDNACollectionsSource,
    GenomeCoordinatesSource as _GenomeCoordinatesSource,
    OligoHybridizationSource as _OligoHybridizationSource,
    PolymeraseExtensionSource as _PolymeraseExtensionSource,
    AnnotationSource as _AnnotationSource,
    AnnotationReport as _AnnotationReport,
    PlannotateAnnotationReport as _PlannotateAnnotationReport,
    ReverseComplementSource as _ReverseComplementSource,
    NCBISequenceSource as _NCBISequenceSource,
    RecombinaseSource as _RecombinaseSource,
    Recombinase as _Recombinase,
)
from Bio.SeqFeature import Location, LocationParserError, SimpleLocation
from Bio.Restriction.Restriction import AbstractCut
import networkx as nx
from typing import List

from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location

from pydna.types import CutSiteType, SubFragmentRepresentationAssembly
from pydna.utils import create_location
from typing import TYPE_CHECKING
import Bio.Restriction as _restr_module
import copy

if TYPE_CHECKING:  # pragma: no cover
    from pydna.dseqrecord import Dseqrecord
    from pydna.primer import Primer
    from pydna.recombinase import RecombinaseCollection, Recombinase


# Thread-local storage for ID strategy
_thread_local = local()


@contextmanager
def id_mode(use_python_internal_id: bool = True):
    """Context manager that is used to determine how ids are assigned to objects when
    mapping them to the OpenCloning data model. If ``use_python_internal_id`` is True,
    the built-in python ``id()`` function is used to assign ids to objects. That function
    produces a unique integer for each object in python, so it's guaranteed to be unique.
    If ``use_python_internal_id`` is False, the object's ``.id`` attribute
    (must be a string integer) is used to assign ids to objects. This is useful
    when the objects already have meaningful ids,
    and you want to keep references to them in ``SourceInput`` objects (which sequences and
    primers are used in a particular source).

    Parameters
    ----------
    use_python_internal_id: bool
        If True, use Python's built-in id() function.
        If False, use the object's .id attribute (must be a string integer).

    Examples
    --------
    >>> from pydna.dseqrecord import Dseqrecord
    >>> from pydna.opencloning_models import get_id, id_mode
    >>> dseqr = Dseqrecord("ATGC")
    >>> dseqr.name = "my_sequence"
    >>> dseqr.id = "123"
    >>> get_id(dseqr) == id(dseqr)
    True
    >>> with id_mode(use_python_internal_id=False):
    ...     get_id(dseqr)
    123
    """
    old_value = getattr(_thread_local, "use_python_internal_id", True)
    _thread_local.use_python_internal_id = use_python_internal_id
    try:
        yield
    finally:
        _thread_local.use_python_internal_id = old_value


def get_id(obj: "Primer" | "Dseqrecord") -> int:
    """Get ID using the current strategy from thread-local storage (see id_mode)
    Parameters
    ----------
    obj: Primer | Dseqrecord
        The object to get the id of

    Returns
    -------
    int: The id of the object

    """
    use_python_internal_id = getattr(_thread_local, "use_python_internal_id", True)
    if use_python_internal_id:
        return id(obj)
    if not isinstance(obj.id, str) or not obj.id.isdigit():
        raise ValueError(
            f"If use_python_internal_id is False, id must be a string representing an integer, "
            f"but object {obj} has an invalid id: {obj.id}"
        )
    return int(obj.id)


class SequenceLocationStr(str):
    """A string representation of a sequence location, genbank-like."""

    @classmethod
    def from_biopython_location(cls, location: Location):
        return cls(format_feature_location(location, None))

    def to_biopython_location(self) -> Location:
        return Location.fromstring(self)

    @classmethod
    def field_validator(cls, v):
        if isinstance(v, str):
            value = cls(v)
            try:
                value.to_biopython_location()
            except LocationParserError as err:
                raise ValueError(f"Location {v!r} is not a valid location") from err
            return value
        raise ValueError(f"Location must be a string or a {cls.__name__}")

    @classmethod
    def __get_pydantic_core_schema__(
        cls,
        source_type,
        handler,
    ) -> core_schema.CoreSchema:
        """Generate Pydantic core schema for SequenceLocationStr."""
        return core_schema.with_info_after_validator_function(
            cls._validate,
            core_schema.str_schema(),
        )

    @classmethod
    def _validate(cls, value: str, info):
        """Validate and create SequenceLocationStr instance."""
        return cls.field_validator(value)

    @classmethod
    def from_start_and_end(
        cls, start: int, end: int, seq_len: int | None = None, strand: int | None = 1
    ):
        return cls.from_biopython_location(create_location(start, end, seq_len, strand))

    def get_ncbi_format_coordinates(self) -> str:
        """Return start, end, strand in the same format as the NCBI eutils API (1-based, inclusive)"""
        return (
            self.to_biopython_location().start + 1,
            self.to_biopython_location().end,
            self.to_biopython_location().strand,
        )


class ConfiguredBaseModel(BaseModel):
    model_config = ConfigDict(
        validate_assignment=True,
        validate_default=True,
        extra="forbid",
        arbitrary_types_allowed=True,
        use_enum_values=True,
        strict=False,
    )
    pass


class TextFileSequence(_TextFileSequence):

    @classmethod
    def from_dseqrecord(cls, dseqr: "Dseqrecord"):
        return cls(
            id=get_id(dseqr),
            sequence_file_format="genbank",
            overhang_crick_3prime=dseqr.seq.ovhg,
            overhang_watson_3prime=dseqr.seq.watson_ovhg,
            file_content=dseqr.format("genbank"),
        )


class PrimerModel(_PrimerModel):

    @classmethod
    def from_primer(cls, primer: "Primer"):
        return cls(
            id=get_id(primer),
            name=primer.name,
            sequence=str(primer.seq),
        )


class SourceInput(ConfiguredBaseModel):
    sequence: object

    @field_validator("sequence")
    @classmethod
    def _validate_sequence_field(cls, value: Any):
        """Separate validation to avoid circular imports."""

        from pydna.dseqrecord import Dseqrecord
        from pydna.primer import Primer

        if isinstance(value, (Dseqrecord, Primer)):
            return value
        module = type(value).__module__
        name = type(value).__name__
        raise TypeError(f"sequence must be Dseqrecord or Primer; got {module}.{name}")

    def to_pydantic_model(self) -> _SourceInput:
        return _SourceInput(sequence=get_id(self.sequence))


class AssemblyFragment(SourceInput):

    left_location: Optional[Location] = Field(default=None)
    right_location: Optional[Location] = Field(default=None)
    reverse_complemented: bool

    @staticmethod
    def from_biopython_location(location: Location | None):
        if location is None:
            return None
        return SequenceLocationStr.from_biopython_location(location)

    def to_pydantic_model(self) -> _AssemblyFragment:
        return _AssemblyFragment(
            sequence=get_id(self.sequence),
            left_location=self.from_biopython_location(self.left_location),
            right_location=self.from_biopython_location(self.right_location),
            reverse_complemented=self.reverse_complemented,
        )


def _source_input_from_model(
    model: _SourceInput | _AssemblyFragment,
    sequences: dict[int, "Dseqrecord"],
    primers: dict[int, "Primer"],
) -> SourceInput | AssemblyFragment:
    """Convert a linkml SourceInput/AssemblyFragment model back to a pydna object."""
    seq_id = model.sequence
    seq_obj = sequences.get(seq_id) or primers.get(seq_id)
    if seq_obj is None:
        raise ValueError(f"Sequence/Primer with id {seq_id} not found")
    if isinstance(model, _AssemblyFragment):
        # It's important to set stranded=False here, because that's
        # how assembly2 creates locations.
        left_loc = (
            Location.fromstring(model.left_location, stranded=False)
            if model.left_location
            else None
        )
        right_loc = (
            Location.fromstring(model.right_location, stranded=False)
            if model.right_location
            else None
        )
        return AssemblyFragment(
            sequence=seq_obj,
            left_location=left_loc,
            right_location=right_loc,
            reverse_complemented=model.reverse_complemented,
        )
    return SourceInput(sequence=seq_obj)


class Source(ConfiguredBaseModel):
    TARGET_MODEL: ClassVar[Type[_Source]] = _Source
    input: list[Union[SourceInput, AssemblyFragment]] = Field(default_factory=list)
    database_id: Optional[int] = None

    @field_serializer("input")
    def serialize_input(
        self, input: list[Union[SourceInput, AssemblyFragment]]
    ) -> list[_SourceInput | _AssemblyFragment]:
        return [fragment.to_pydantic_model() for fragment in input]

    def to_pydantic_model(self, seq_id: int):
        model_dict = self.model_dump()
        model_dict["id"] = seq_id
        return self.TARGET_MODEL(**model_dict)

    @classmethod
    def _get_deserialization_overrides(cls, model, sequences, primers) -> dict:
        """Return field overrides for from_pydantic_model. Subclasses override this
        to handle fields that require custom deserialization logic."""
        return {}

    @classmethod
    def from_pydantic_model(
        cls,
        model: _Source,
        sequences: dict[int, "Dseqrecord"],
        primers: dict[int, "Primer"],
    ) -> "Source":
        """Convert an opencloning_linkml Source model back to a pydna Source."""
        overrides = cls._get_deserialization_overrides(model, sequences, primers)
        kwargs = {}
        for field_name in cls.__pydantic_fields__:
            if field_name in overrides:
                kwargs[field_name] = overrides[field_name]
            elif field_name == "input":
                kwargs["input"] = [
                    _source_input_from_model(inp, sequences, primers)
                    for inp in (model.input or [])
                ]
            else:
                kwargs[field_name] = getattr(model, field_name)
        return cls(**kwargs)

    def to_unserialized_dict(self):
        """
        Converts into a dictionary without serializing the fields.
        This is used to be able to recast.
        """
        return {field: getattr(self, field) for field in self.__pydantic_fields__}

    def add_to_history_graph(self, history_graph: nx.DiGraph, seq: "Dseqrecord"):
        """
        Add the source to the history graph.

        It does not use the get_id function, because it just uses it to have unique identifiers
        for graph nodes, not to store them anywhere.
        """
        from pydna.dseqrecord import Dseqrecord

        history_graph.add_node(id(seq), label=f"{seq.name} ({repr(seq)})")
        history_graph.add_node(id(self), label=str(self.TARGET_MODEL.__name__))
        history_graph.add_edge(id(seq), id(self))
        for fragment in self.input:
            fragment_seq = fragment.sequence
            # This could be a Primer as well, which doesn't have a source
            if isinstance(fragment_seq, Dseqrecord) and fragment_seq.source is not None:
                fragment_seq.source.add_to_history_graph(history_graph, fragment_seq)
            else:
                history_graph.add_node(
                    id(fragment_seq),
                    label=f"{fragment_seq.name} ({repr(fragment_seq)})",
                )
            history_graph.add_edge(id(self), id(fragment_seq))

    def history_string(self, seq: "Dseqrecord"):
        """
        Returns a string representation of the cloning history of the sequence.
        See dseqrecord.history() for examples.
        """
        history_graph = nx.DiGraph()
        self.add_to_history_graph(history_graph, seq)
        return "\n".join(
            nx.generate_network_text(history_graph, with_labels=True, sources=[id(seq)])
        )

    def _replay_products(self, **kwargs) -> list["Dseqrecord"]:
        """Replay the cloning operation and return all possible products.

        Subclasses must override this to provide the replay logic.
        """
        raise NotImplementedError(
            f"_replay_products() not implemented for {type(self).__name__}"
        )

    def _validate_result_in_products(
        self, result: "Dseqrecord", products: list["Dseqrecord"]
    ) -> None:
        if not any((p.seq == result.seq) and (p.source == self) for p in products):
            raise ValueError(
                f"Result sequence does not match any of the {len(products)} "
                f"product(s) from {type(self).__name__}"
            )

    def _find_product_by_seguid(
        self, result: "Dseqrecord", products: list["Dseqrecord"]
    ) -> "Dseqrecord":
        """Find and return the product whose seguid matches result's seguid."""
        target_seguid = result.seq.seguid()
        for p in products:
            if p.seq.seguid() == target_seguid:
                return p
        raise ValueError(
            f"No product with matching seguid among {len(products)} "
            f"product(s) from {type(self).__name__}"
        )

    def validate(self, result: "Dseqrecord") -> None:
        """Replay the cloning operation and verify the result sequence matches.

        Raises ValueError if the result doesn't match any expected product.
        Raises NotImplementedError for source types that don't support validation yet.
        Sources with no inputs (external imports) return silently.
        """
        if not self.input:
            return
        self._validate_result_in_products(result, self._replay_products())

    def normalize(self, result: "Dseqrecord") -> "Dseqrecord":
        """Replay the cloning operation and return the product matching result by seguid.

        Uses seguid comparison which is invariant to rotation and orientation,
        so it finds the correct product even if inputs were rotated or reverse-complemented.
        Copies name and id from the original result onto the returned product.
        """
        if not self.input:
            return result
        product = self._find_product_by_seguid(result, self._replay_products())
        product.name = result.name
        product.id = result.id
        return product


class AssemblySource(Source):
    circular: bool
    input: list[Union[SourceInput, AssemblyFragment]] = Field(
        default_factory=list, min_length=1
    )

    TARGET_MODEL: ClassVar[Type[_AssemblySource]] = _AssemblySource

    @classmethod
    def from_subfragment_representation(
        cls,
        assembly: SubFragmentRepresentationAssembly,
        fragments: list["Dseqrecord"],
        is_circular: bool,
    ):

        input_list = []
        for f_index, loc1, loc2 in assembly:
            input_list.append(
                AssemblyFragment(
                    sequence=fragments[abs(f_index) - 1],
                    left_location=loc1,
                    right_location=loc2,
                    reverse_complemented=f_index < 0,
                )
            )

        return AssemblySource(input=input_list, circular=is_circular)

    def _get_input_sequences(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        """Return Dseqrecord inputs (excludes Primers), preserving order, and handling insertion assemblies (removing the last fragment if the first and last are the same)."""
        from pydna.dseqrecord import Dseqrecord

        seqs: list[Dseqrecord] = []
        for inp in self.input:
            if isinstance(inp.sequence, Dseqrecord):
                seqs.append(inp.sequence)

        # This is necessary for insertion assemblies, where the first and last fragment are the same.
        if handle_insertion and len(seqs) > 1 and seqs[0] == seqs[-1]:
            return seqs[:-1]
        return seqs

    def _get_input_primers(self) -> list["Primer"]:
        """Return Primer inputs, preserving order."""
        from pydna.primer import Primer

        return [inp.sequence for inp in self.input if isinstance(inp.sequence, Primer)]

    def _minimal_assembly_overlap(self) -> int:
        """Return the smallest overlap length across all AssemblyFragment locations."""
        overlaps: list[int] = []
        for inp in self.input:
            if not isinstance(inp, AssemblyFragment):
                continue
            if inp.left_location is not None:
                overlaps.append(len(inp.left_location))
            if inp.right_location is not None:
                overlaps.append(len(inp.right_location))
        if not overlaps:
            raise ValueError("Assembly is not complete")
        return min(overlaps)

    def validate(self, result: "Dseqrecord") -> None:
        products = self._replay_products(True)
        try:
            self._validate_result_in_products(result, products)
        except ValueError:
            products = self._replay_products(False)
            self._validate_result_in_products(result, products)

    def normalize(self, result: "Dseqrecord") -> "Dseqrecord":
        products = self._replay_products(True)
        try:
            product = self._find_product_by_seguid(result, products)
        except ValueError:
            products = self._replay_products(False)
            product = self._find_product_by_seguid(result, products)
        product.name = result.name
        product.id = result.id
        return product

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        super()._replay_products(handle_insertion=handle_insertion)


_empty_input_list = Field(
    default_factory=list,
    min_length=0,
    max_length=0,
    description="Must always be an empty list.",
)


class DatabaseSource(Source):
    TARGET_MODEL: ClassVar[Type[_DatabaseSource]] = _DatabaseSource
    input: list[Union[SourceInput, AssemblyFragment]] = _empty_input_list

    database_id: int


class UploadedFileSource(Source):
    TARGET_MODEL: ClassVar[Type[_UploadedFileSource]] = _UploadedFileSource
    input: list[Union[SourceInput, AssemblyFragment]] = _empty_input_list

    file_name: str
    index_in_file: int
    sequence_file_format: str


class RepositoryIdSource(Source):
    TARGET_MODEL: ClassVar[Type[_RepositoryIdSource]] = _RepositoryIdSource
    input: list[Union[SourceInput, AssemblyFragment]] = _empty_input_list

    repository_id: str
    # location: Location


class RepositoryIdSourceWithSequenceFileUrl(RepositoryIdSource):
    """
    Auxiliary class to avoid code duplication in the sources that have
    a sequence file url.
    """

    input: list[Union[SourceInput, AssemblyFragment]] = _empty_input_list
    sequence_file_url: Optional[str] = None


class AddgeneIdSource(RepositoryIdSourceWithSequenceFileUrl):
    TARGET_MODEL: ClassVar[Type[_AddgeneIdSource]] = _AddgeneIdSource
    addgene_sequence_type: Optional[AddgeneSequenceType] = None


class BenchlingUrlSource(RepositoryIdSource):
    TARGET_MODEL: ClassVar[Type[_BenchlingUrlSource]] = _BenchlingUrlSource


class SnapGenePlasmidSource(RepositoryIdSource):
    TARGET_MODEL: ClassVar[Type[_SnapGenePlasmidSource]] = _SnapGenePlasmidSource


class EuroscarfSource(RepositoryIdSource):
    TARGET_MODEL: ClassVar[Type[_EuroscarfSource]] = _EuroscarfSource


class WekWikGeneIdSource(RepositoryIdSourceWithSequenceFileUrl):
    TARGET_MODEL: ClassVar[Type[_WekWikGeneIdSource]] = _WekWikGeneIdSource


class SEVASource(RepositoryIdSourceWithSequenceFileUrl):
    TARGET_MODEL: ClassVar[Type[_SEVASource]] = _SEVASource


class IGEMSource(RepositoryIdSourceWithSequenceFileUrl):
    TARGET_MODEL: ClassVar[Type[_IGEMSource]] = _IGEMSource


class OpenDNACollectionsSource(RepositoryIdSourceWithSequenceFileUrl):
    TARGET_MODEL: ClassVar[Type[_OpenDNACollectionsSource]] = _OpenDNACollectionsSource


class NCBISequenceSource(RepositoryIdSource):
    TARGET_MODEL: ClassVar[Type[_NCBISequenceSource]] = _NCBISequenceSource
    coordinates: SimpleLocation | None = None

    @classmethod
    def _get_deserialization_overrides(
        cls, model: "_NCBISequenceSource", sequences, primers
    ):
        return {
            "coordinates": (
                Location.fromstring(model.coordinates) if model.coordinates else None
            )
        }


class GenomeCoordinatesSource(NCBISequenceSource):
    TARGET_MODEL: ClassVar[Type[_GenomeCoordinatesSource]] = _GenomeCoordinatesSource

    assembly_accession: Optional[str] = None
    locus_tag: Optional[str] = None
    gene_id: Optional[int] = None
    coordinates: SimpleLocation

    @field_serializer("coordinates")
    def serialize_coordinates(self, coordinates: SimpleLocation) -> str:
        return SequenceLocationStr.from_biopython_location(coordinates)


class RestrictionAndLigationSource(AssemblySource):
    restriction_enzymes: list[AbstractCut]

    TARGET_MODEL: ClassVar[Type[_RestrictionAndLigationSource]] = (
        _RestrictionAndLigationSource
    )

    @field_serializer("restriction_enzymes")
    def serialize_restriction_enzymes(
        self, restriction_enzymes: list[AbstractCut]
    ) -> list[str]:
        return [str(enzyme) for enzyme in restriction_enzymes]

    @classmethod
    def _get_deserialization_overrides(cls, model, sequences, primers):
        return {
            "restriction_enzymes": [
                getattr(_restr_module, name) for name in model.restriction_enzymes
            ]
        }

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import restriction_ligation_assembly

        return restriction_ligation_assembly(
            self._get_input_sequences(handle_insertion), self.restriction_enzymes
        )


class GibsonAssemblySource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_GibsonAssemblySource]] = _GibsonAssemblySource

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import gibson_assembly

        return gibson_assembly(
            self._get_input_sequences(handle_insertion),
            limit=self._minimal_assembly_overlap(),
        )


class InFusionSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_InFusionSource]] = _InFusionSource

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import in_fusion_assembly

        return in_fusion_assembly(
            self._get_input_sequences(handle_insertion),
            limit=self._minimal_assembly_overlap(),
        )


class OverlapExtensionPCRLigationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_OverlapExtensionPCRLigationSource]] = (
        _OverlapExtensionPCRLigationSource
    )

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import fusion_pcr_assembly

        return fusion_pcr_assembly(
            self._get_input_sequences(handle_insertion),
            limit=self._minimal_assembly_overlap(),
        )


class InVivoAssemblySource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_InVivoAssemblySource]] = _InVivoAssemblySource

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import in_vivo_assembly

        return in_vivo_assembly(
            self._get_input_sequences(handle_insertion),
            limit=self._minimal_assembly_overlap(),
        )


class LigationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_LigationSource]] = _LigationSource

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import ligation_assembly

        overlap = self._minimal_assembly_overlap()
        return ligation_assembly(
            self._get_input_sequences(handle_insertion),
            allow_blunt=overlap == 0,
            allow_partial_overlap=True,
        )


class GatewaySource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_GatewaySource]] = _GatewaySource
    reaction_type: GatewayReactionType
    greedy: bool = Field(default=False)

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import gateway_assembly

        return gateway_assembly(
            self._get_input_sequences(handle_insertion),
            self.reaction_type,
            self.greedy,
        )


class HomologousRecombinationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_HomologousRecombinationSource]] = (
        _HomologousRecombinationSource
    )

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import (
            homologous_recombination_integration,
            homologous_recombination_excision,
        )

        seqs = self._get_input_sequences(handle_insertion)
        limit = self._minimal_assembly_overlap()
        if len(seqs) == 1:
            return homologous_recombination_excision(seqs[0], limit)
        else:
            return homologous_recombination_integration(seqs[0], seqs[1:], limit)


class CRISPRSource(HomologousRecombinationSource):
    TARGET_MODEL: ClassVar[Type[_CRISPRSource]] = _CRISPRSource

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import crispr_integration

        seqs = self._get_input_sequences(handle_insertion)
        guides = self._get_input_primers()
        limit = self._minimal_assembly_overlap()
        return crispr_integration(seqs[0], seqs[1:], guides, limit)


class CreLoxRecombinationSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_CreLoxRecombinationSource]] = (
        _CreLoxRecombinationSource
    )

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import cre_lox_integration, cre_lox_excision

        seqs = self._get_input_sequences(handle_insertion)
        if len(seqs) == 1:
            return cre_lox_excision(seqs[0])
        else:
            return cre_lox_integration(seqs[0], seqs[1:])


class PCRSource(AssemblySource):
    TARGET_MODEL: ClassVar[Type[_PCRSource]] = _PCRSource
    add_primer_features: bool = Field(default=False)

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import pcr_assembly

        seqs = self._get_input_sequences(handle_insertion)
        primers = self._get_input_primers()
        limit = self._minimal_assembly_overlap()
        return pcr_assembly(
            seqs[0],
            primers[0],
            primers[1],
            limit=limit,
            add_primer_features=self.add_primer_features,
        )


class RecombinaseSource(AssemblySource):

    TARGET_MODEL: ClassVar[Type[_RecombinaseSource]] = _RecombinaseSource
    recombinases: object

    @field_validator("recombinases")
    @classmethod
    def _validate_recombinases(cls, value: "RecombinaseCollection" | "Recombinase"):
        from pydna.recombinase import RecombinaseCollection, Recombinase

        if isinstance(value, RecombinaseCollection) or isinstance(value, Recombinase):
            return value
        raise ValueError(
            "recombinases must be a RecombinaseCollection or a Recombinase"
        )

    @field_serializer("recombinases")
    def serialize_recombinases(
        self, recombinases: "RecombinaseCollection" | "Recombinase"
    ) -> list[_Recombinase]:
        return recombinases.to_opencloning_model()

    @classmethod
    def _get_deserialization_overrides(
        cls, model: "_RecombinaseSource", sequences, primers
    ):
        from pydna.recombinase import (
            RecombinaseCollection,
            Recombinase as _RecombinaseClass,
        )

        recomb_objs = [
            _RecombinaseClass(
                site1=r.site1,
                site2=r.site2,
                site1_name=r.site1_name,
                site2_name=r.site2_name,
                name=r.name,
            )
            for r in model.recombinases
        ]
        return {
            "recombinases": (
                recomb_objs[0]
                if len(recomb_objs) == 1
                else RecombinaseCollection(recomb_objs)
            )
        }

    def _replay_products(self, handle_insertion: bool = True) -> list["Dseqrecord"]:
        from pydna.assembly2 import recombinase_integration, recombinase_excision

        seqs = self._get_input_sequences(handle_insertion)
        if len(seqs) == 1:
            return recombinase_excision(seqs[0], self.recombinases)
        else:
            return recombinase_integration(seqs[0], seqs[1:], self.recombinases)


class SequenceCutSource(Source):
    left_edge: CutSiteType | None
    right_edge: CutSiteType | None
    input: list[Union[SourceInput, AssemblyFragment]] = Field(
        default_factory=list, min_length=1, max_length=1
    )

    @property
    def TARGET_MODEL(self):
        return (
            _RestrictionEnzymeDigestionSource
            if self._has_enzyme()
            else _SequenceCutSource
        )

    @field_serializer("left_edge", "right_edge")
    def serialize_cut_site(
        self, cut_site: CutSiteType | None
    ) -> _RestrictionSequenceCut | _SequenceCut | None:
        return self._cutsite_to_model(cut_site)

    @staticmethod
    def _cutsite_to_model(cut_site: CutSiteType | None):
        if cut_site is None:
            return None
        watson, overhang = cut_site[0]
        enzyme_or_none = cut_site[1]
        if isinstance(enzyme_or_none, AbstractCut):
            return _RestrictionSequenceCut(
                cut_watson=watson,
                overhang=overhang,
                restriction_enzyme=str(enzyme_or_none),
            )
        return _SequenceCut(cut_watson=watson, overhang=overhang)

    @classmethod
    def from_parent(
        cls, parent: "Dseqrecord", left_edge: CutSiteType, right_edge: CutSiteType
    ):
        return cls(
            input=[SourceInput(sequence=parent)],
            left_edge=left_edge,
            right_edge=right_edge,
        )

    def _has_enzyme(self) -> bool:
        def has_enzyme(edge):
            return edge is not None and isinstance(edge[1], AbstractCut)

        return has_enzyme(self.left_edge) or has_enzyme(self.right_edge)

    def _replay_products(self) -> list["Dseqrecord"]:
        parent = self.input[0].sequence
        if self._has_enzyme():
            enzymes = set()
            for edge in (self.left_edge, self.right_edge):
                if edge is not None and isinstance(edge[1], AbstractCut):
                    enzymes.add(edge[1])
            return list(parent.cut(*enzymes))
        return [parent.apply_cut(self.left_edge, self.right_edge)]

    def normalize(self, result: "Dseqrecord") -> "Dseqrecord":
        parent = self.input[0].sequence
        if self._has_enzyme() or parent.seq == result.seq:
            return super().normalize(result)
        raise ValueError(
            "SequenceCutSource must have an enzyme to be normalized if sequences differ"
        )

    @classmethod
    def _get_deserialization_overrides(
        cls, model: "_RestrictionEnzymeDigestionSource", sequences, primers
    ):
        def _model_to_cutsite(edge_model):
            if edge_model is None:
                return None
            enzyme = None
            if isinstance(edge_model, _RestrictionSequenceCut):
                enzyme = getattr(_restr_module, edge_model.restriction_enzyme)
            return ((edge_model.cut_watson, edge_model.overhang), enzyme)

        return {
            "left_edge": _model_to_cutsite(model.left_edge),
            "right_edge": _model_to_cutsite(model.right_edge),
        }


class OligoHybridizationSource(Source):
    TARGET_MODEL: ClassVar[Type[_OligoHybridizationSource]] = _OligoHybridizationSource
    input: list[Union[SourceInput]] = Field(
        default_factory=list, min_length=2, max_length=2
    )
    overhang_crick_3prime: Optional[int] = None

    def _replay_products(self) -> list["Dseqrecord"]:
        from pydna.oligonucleotide_hybridization import oligonucleotide_hybridization

        fwd = self.input[0].sequence
        rvs = self.input[1].sequence
        annealing = min(len(fwd.seq), len(rvs.seq)) - abs(
            self.overhang_crick_3prime or 0
        )
        return oligonucleotide_hybridization(fwd, rvs, annealing)


class PolymeraseExtensionSource(Source):
    TARGET_MODEL: ClassVar[Type[_PolymeraseExtensionSource]] = (
        _PolymeraseExtensionSource
    )
    input: list[Union[SourceInput]] = Field(
        default_factory=list, min_length=1, max_length=1
    )

    def _replay_products(self) -> list["Dseqrecord"]:
        from pydna.dseqrecord import Dseqrecord

        inp = self.input[0].sequence
        if not isinstance(inp, Dseqrecord):
            raise ValueError(
                f"PolymeraseExtensionSource input must be a Dseqrecord, got {type(inp)}"
            )
        result = copy.deepcopy(inp)
        result.seq = result.seq.fill_in()
        result.source = self
        return [result]


class AnnotationSource(Source):
    TARGET_MODEL: ClassVar[Type[_AnnotationSource]] = _AnnotationSource
    input: list[Union[SourceInput]] = Field(
        default_factory=list, min_length=1, max_length=1
    )
    annotation_tool: AnnotationTool
    annotation_tool_version: Optional[str] = None
    annotation_report: Optional[
        list[_AnnotationReport | _PlannotateAnnotationReport]
    ] = None

    def validate(self, result: "Dseqrecord") -> None:
        """Just validates that there is a single input, and that its sequence is the same as the result."""
        from pydna.dseqrecord import Dseqrecord

        if not isinstance(self.input[0].sequence, Dseqrecord):
            raise ValueError(
                f"AnnotationSource input must be a Dseqrecord, got {type(self.input[0].sequence)}"
            )
        if self.input[0].sequence.seq != result.seq:
            raise ValueError("AnnotationSource input sequence does not match result")

    def normalize(self, result: "Dseqrecord") -> "Dseqrecord":
        from pydna.dseqrecord import Dseqrecord

        input_sequence = self.input[0].sequence
        if not isinstance(input_sequence, Dseqrecord):
            raise ValueError(
                f"AnnotationSource input must be a Dseqrecord, got {type(input_sequence)}"
            )
        if input_sequence.seq.seguid() != result.seq.seguid():
            raise ValueError("AnnotationSource input sequence does not match result")
        if input_sequence.seq == result.seq:
            return copy.deepcopy(result)
        if input_sequence.circular:
            # synced already handles reverse complementing
            out_seq = result.synced(input_sequence.seq)
        else:
            out_seq = result.reverse_complement()
        out_source = copy.deepcopy(self)
        # Here we remove the annotation report, since it referenced original coordinates
        out_source.annotation_report = []
        out_seq.source = out_source
        return out_seq


class ReverseComplementSource(Source):
    TARGET_MODEL: ClassVar[Type[_ReverseComplementSource]] = _ReverseComplementSource
    input: list[Union[SourceInput]] = Field(
        default_factory=list, min_length=1, max_length=1
    )

    def _replay_products(self) -> list["Dseqrecord"]:
        out_seq = self.input[0].sequence.reverse_complement()
        out_seq.source = self
        return [out_seq]


def read_dseqrecord_from_text_file_sequence(seq: TextFileSequence) -> Dseqrecord:
    from pydna.parsers import parse
    from pydna.dseq import Dseq

    out_dseq_record = parse(seq.file_content, ds=True)[0]
    if seq.overhang_watson_3prime != 0 or seq.overhang_crick_3prime != 0:
        out_dseq_record.seq = Dseq.from_full_sequence_and_overhangs(
            str(out_dseq_record.seq),
            seq.overhang_crick_3prime,
            seq.overhang_watson_3prime,
        )
    out_dseq_record.id = str(seq.id)
    return out_dseq_record


class CloningStrategy(_BaseCloningStrategy):

    # For now, we don't add anything, but the classes will not have the new
    # methods if this is used
    # It will be used for validation for now
    primers: Optional[List[PrimerModel]] = Field(
        default_factory=list,
        description="""The primers that are used in the cloning strategy""",
        json_schema_extra={
            "linkml_meta": {"alias": "primers", "domain_of": ["CloningStrategy"]}
        },
    )

    def add_primer(self, primer: "Primer"):
        existing_ids = {seq.id for seq in self.primers}
        if get_id(primer) in existing_ids:
            return
        self.primers.append(PrimerModel.from_primer(primer))

    def add_dseqrecord(self, dseqr: "Dseqrecord"):
        from pydna.dseqrecord import Dseqrecord

        existing_ids = {seq.id for seq in self.sequences}
        if get_id(dseqr) in existing_ids:
            return
        self.sequences.append(TextFileSequence.from_dseqrecord(dseqr))
        if dseqr.source is not None:
            self.sources.append(dseqr.source.to_pydantic_model(get_id(dseqr)))
            this_source: Source = dseqr.source
            for source_input in this_source.input:
                if isinstance(source_input.sequence, Dseqrecord):
                    self.add_dseqrecord(source_input.sequence)
                else:
                    self.add_primer(source_input.sequence)
        else:
            self.sources.append(_ManuallyTypedSource(id=get_id(dseqr), input=[]))

    def reassign_ids(self):
        all_ids = (
            {seq.id for seq in self.sequences}
            | {source.id for source in self.sources}
            | {primer.id for primer in self.primers}
        )
        id_mappings = {id: i + 1 for i, id in enumerate(sorted(all_ids))}
        for seq in self.sequences:
            seq.id = id_mappings[seq.id]
        for primer in self.primers:
            primer.id = id_mappings[primer.id]
        for source in self.sources:
            source.id = id_mappings[source.id]
            for assembly_fragment in source.input:
                assembly_fragment.sequence = id_mappings[assembly_fragment.sequence]

    @classmethod
    def from_dseqrecords(cls, dseqrs: list["Dseqrecord"], description: str = ""):
        cloning_strategy = cls(sources=[], sequences=[], description=description)
        for dseqr in dseqrs:
            cloning_strategy.add_dseqrecord(dseqr)
        return cloning_strategy

    def model_dump_json(self, *args, **kwargs):
        if getattr(_thread_local, "use_python_internal_id", True):
            # Make a deep copy of the cloning strategy and reassign ids
            cs = self.__deepcopy__()
            cs.reassign_ids()
            return super(CloningStrategy, cs).model_dump_json(*args, **kwargs)
        return super().model_dump_json(*args, **kwargs)

    def model_dump(self, *args, **kwargs):
        if getattr(_thread_local, "use_python_internal_id", True):
            cs = self.__deepcopy__()
            cs.reassign_ids()
            return super(CloningStrategy, cs).model_dump(*args, **kwargs)
        return super().model_dump(*args, **kwargs)

    def get_ids_of_sequences_that_are_not_inputs(self) -> list[int]:
        """
        Get the ids of the sequences that are not inputs to any source.
        """
        return [
            seq.id
            for seq in self.sequences
            if seq.id
            not in [input.sequence for source in self.sources for input in source.input]
        ]

    def get_ids_of_sequences_that_are_inputs(self) -> list[int]:
        """
        Get the ids of the sequences that are inputs to any source.
        """
        return [
            seq.id
            for seq in self.sequences
            if seq.id
            in [input.sequence for source in self.sources for input in source.input]
        ]

    def to_dseqrecords(self) -> list["Dseqrecord"]:
        """
        Convert the cloning strategy to a list of Dseqrecord objects.

        Walks the dependency graph from terminal sequences (those not used as inputs
        to any source) upward, reconstructing pydna Source objects and attaching them
        to the parsed Dseqrecord objects.
        """
        from pydna.primer import Primer as _PrimerCls

        # Build primer lookup
        primer_by_id: dict[int, _PrimerCls] = {}
        for p in self.primers or []:
            primer_by_id[p.id] = _PrimerCls(p.sequence, name=p.name)

        # Build sequence lookup by parsing TextFileSequence genbank content
        seq_by_id: dict[int, Dseqrecord] = {}
        for text_seq in self.sequences:
            if not isinstance(text_seq, _TextFileSequence):
                raise NotImplementedError(
                    f"Only works with TextFileSequence objects, got {type(text_seq)}"
                )
            seq_by_id[text_seq.id] = read_dseqrecord_from_text_file_sequence(text_seq)

        # Build source lookup (source.id == output sequence id)
        source_by_id: dict[int, _Source] = {s.id: s for s in self.sources}

        # Resolve sources recursively
        resolved: set[int] = set()

        def resolve(seq_id: int):
            if seq_id in resolved:
                return
            resolved.add(seq_id)

            source_model = source_by_id.get(seq_id)
            if source_model is None:
                raise ValueError(f"Missing source for sequence {seq_id}")

            # Recursively resolve input sequences first
            for inp in source_model.input or []:
                if inp.sequence in seq_by_id:
                    resolve(inp.sequence)

            dseqr = seq_by_id[seq_id]

            if isinstance(source_model, _ManuallyTypedSource):
                dseqr.source = None
                return

            pydna_cls = _TARGET_MODEL_REGISTRY.get(type(source_model))
            if pydna_cls is None:  # pragma: no cover
                raise ValueError(f"Unknown source model type: {type(source_model)}")

            dseqr.source = pydna_cls.from_pydantic_model(
                source_model, seq_by_id, primer_by_id
            )

        terminal_ids = self.get_ids_of_sequences_that_are_not_inputs()
        for sid in terminal_ids:
            resolve(sid)

        return [seq_by_id[sid] for sid in terminal_ids]

    def normalize(self) -> "CloningStrategy":
        """
        Normalize the cloning strategy by resolving the sources and sequences.
        """
        newseqrs = [s.normalize_history() for s in self.to_dseqrecords()]
        return self.__class__.from_dseqrecords(newseqrs, self.description)

    def validate(self) -> None:
        """
        Validate the cloning strategy by resolving the sources and sequences.
        """
        for seq in self.to_dseqrecords():
            seq.validate_history()


def _build_target_model_registry() -> dict[type, type]:
    """Build a mapping from opencloning_linkml Source types to pydna Source types."""
    registry = {}

    def _register(cls):
        target = cls.__dict__.get("TARGET_MODEL")
        if target is not None and not isinstance(target, property):
            registry[target] = cls
        for sub in cls.__subclasses__():
            _register(sub)

    _register(Source)
    # SequenceCutSource uses @property for TARGET_MODEL, register manually
    registry[_SequenceCutSource] = SequenceCutSource
    registry[_RestrictionEnzymeDigestionSource] = SequenceCutSource
    return registry


_TARGET_MODEL_REGISTRY: dict[type, type] = _build_target_model_registry()
