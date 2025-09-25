# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Optional, Union
from pydantic_core import core_schema

from pydantic import BaseModel, ConfigDict, Field

from opencloning_linkml.datamodel import (
    CloningStrategy as _BaseCloningStrategy,
    Primer as PrimerModel,
    Source as _Source,
    TextFileSequence as _TextFileSequence,
    AssemblySource as _AssemblySource,
    SourceInput as _SourceInput,
    AssemblyFragment as _AssemblyFragment,
    ManuallyTypedSource as _ManuallyTypedSource,
    RestrictionAndLigationSource as _RestrictionAndLigationSource,
    GibsonAssemblySource as _GibsonAssemblySource,
    RestrictionEnzymeDigestionSource as _RestrictionEnzymeDigestionSource,
    RestrictionSequenceCut as _RestrictionSequenceCut,
)
from Bio.SeqFeature import Location, LocationParserError
from Bio.Restriction.Restriction import AbstractCut
from typing import List

from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location

from pydna.types import CutSiteType
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pydna.dseqrecord import Dseqrecord


class SequenceLocationStr(str):
    """A string representation of a sequence location, genbank-like."""

    # TODO: this should handle origin-spanning simple locations (splitted)
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
            except LocationParserError:
                raise ValueError(f'Location "{v}" is not a valid location')
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
            id=id(dseqr),
            sequence_file_format="genbank",
            overhang_crick_3prime=dseqr.seq.ovhg,
            overhang_watson_3prime=dseqr.seq.watson_ovhg(),
            file_content=dseqr.format("genbank"),
        )


class SourceInput(ConfiguredBaseModel):
    sequence: "Dseqrecord"

    def to_pydantic_model(self) -> _SourceInput:
        return _SourceInput(sequence=id(self.sequence))


class AssemblyFragment(SourceInput):
    """
    Represents a fragment in an assembly
    """

    left_location: Optional[Location] = Field(default=None)
    right_location: Optional[Location] = Field(default=None)
    reverse_complemented: bool = Field(
        default=...,
    )

    def to_pydantic_model(self) -> _AssemblyFragment:
        return _AssemblyFragment(
            sequence=id(self.sequence),
            left_location=SequenceLocationStr.from_biopython_location(
                self.left_location
            ),
            right_location=SequenceLocationStr.from_biopython_location(
                self.right_location
            ),
            reverse_complemented=self.reverse_complemented,
        )


class Source(ConfiguredBaseModel):
    input: Optional[list[Union[SourceInput, AssemblyFragment]]] = Field(default=None)


class AssemblySource(Source):
    circular: bool

    def to_pydantic_model(self, seq_id) -> _AssemblySource:

        return _AssemblySource(
            id=seq_id,
            input=[fragment.to_pydantic_model() for fragment in self.input],
            circular=self.circular,
        )


class RestrictionAndLigationSource(AssemblySource):
    restriction_enzymes: list[AbstractCut]

    def to_pydantic_model(self, seq_id) -> _RestrictionAndLigationSource:
        parent_model = super().to_pydantic_model(seq_id)
        return _RestrictionAndLigationSource(
            **{k: v for k, v in parent_model.model_dump().items() if k != "type"},
            restriction_enzymes=[str(enzyme) for enzyme in self.restriction_enzymes],
        )


class GibsonAssemblySource(AssemblySource):
    def to_pydantic_model(self, seq_id) -> _GibsonAssemblySource:
        parent_model = super().to_pydantic_model(seq_id)
        return _GibsonAssemblySource(
            **{k: v for k, v in parent_model.model_dump().items() if k != "type"},
        )


class RestrictionEnzymeDigestionSource(Source):
    left_edge: CutSiteType | None
    right_edge: CutSiteType | None

    @classmethod
    def from_parent(
        cls, parent: "Dseqrecord", left_edge: CutSiteType, right_edge: CutSiteType
    ):
        return cls(
            input=[SourceInput(sequence=parent)],
            left_edge=left_edge,
            right_edge=right_edge,
        )

    def to_pydantic_model(self, seq_id) -> _RestrictionEnzymeDigestionSource:
        left_edge = (
            None
            if self.left_edge is None
            else _RestrictionSequenceCut(
                cut_watson=self.left_edge[0][0],
                overhang=self.left_edge[0][1],
                restriction_enzyme=str(self.left_edge[1]),
            )
        )
        right_edge = (
            None
            if self.right_edge is None
            else _RestrictionSequenceCut(
                cut_watson=self.right_edge[0][0],
                overhang=self.right_edge[0][1],
                restriction_enzyme=str(self.right_edge[1]),
            )
        )
        return _RestrictionEnzymeDigestionSource(
            id=seq_id,
            input=[fragment.to_pydantic_model() for fragment in self.input],
            left_edge=left_edge,
            right_edge=right_edge,
        )


class CloningStrategy(_BaseCloningStrategy):

    # For now, we don't add anything, but the classes will not have the new methods if this is used
    # It will be used for validation for now
    primers: Optional[List[PrimerModel]] = Field(
        default_factory=list,
        description="""The primers that are used in the cloning strategy""",
        json_schema_extra={
            "linkml_meta": {"alias": "primers", "domain_of": ["CloningStrategy"]}
        },
    )

    def add_primer(self, primer: PrimerModel):
        if primer in self.primers:
            return
        primer.id = self.next_id()
        self.primers.append(primer)

    def next_id(self):
        return (
            max([s.id for s in self.sources + self.sequences + self.primers], default=0)
            + 1
        )

    def add_source_and_sequence(self, source: _Source, sequence: TextFileSequence):
        if source in self.sources:
            if sequence not in self.sequences:
                raise ValueError(
                    f"Source {source.id} already exists in the cloning strategy, but sequence {sequence.id} it's not its output."
                )
            return
        new_id = self.next_id()
        source.id = new_id
        self.sources.append(source)
        sequence.id = new_id
        self.sequences.append(sequence)

    def all_children_source_ids(
        self, source_id: int, source_children: list | None = None
    ) -> list[int]:
        """Returns the ids of all source children ids of a source"""
        source = next(s for s in self.sources if s.id == source_id)
        if source_children is None:
            source_children = []

        sources_that_take_output_as_input = [
            s for s in self.sources if source.id in [inp.sequence for inp in s.input]
        ]
        new_source_ids = [s.id for s in sources_that_take_output_as_input]

        source_children.extend(new_source_ids)
        for new_source_id in new_source_ids:
            self.all_children_source_ids(new_source_id, source_children)
        return source_children

    def add_dseqrecord(self, dseqr: "Dseqrecord"):
        existing_ids = {seq.id for seq in self.sequences}
        if id(dseqr) in existing_ids:
            return
        self.sequences.append(TextFileSequence.from_dseqrecord(dseqr))
        if dseqr.source is not None:
            self.sources.append(dseqr.source.to_pydantic_model(id(dseqr)))
            this_source: Source = dseqr.source
            for source_input in this_source.input:
                self.add_dseqrecord(source_input.sequence)
        else:
            self.sources.append(
                _ManuallyTypedSource(id=id(dseqr), input=[], user_input="A")
            )

    def reassign_ids(self):
        all_ids = (
            {seq.id for seq in self.sequences}
            | {source.id for source in self.sources}
            | {primer.id for primer in self.primers}
        )
        id_mappings = {id: i + 1 for i, id in enumerate(sorted(all_ids))}
        for seq in self.sequences:
            seq.id = id_mappings[seq.id]
        for source in self.sources:
            source.id = id_mappings[source.id]
            for assembly_fragment in source.input:
                assembly_fragment.sequence = id_mappings[assembly_fragment.sequence]
        for primer in self.primers:
            primer.id = id_mappings[primer.id]

    @classmethod
    def from_dseqrecord(cls, dseqr: "Dseqrecord", description: str = ""):
        cloning_strategy = cls(sources=[], sequences=[], description=description)
        cloning_strategy.add_dseqrecord(dseqr)
        return cloning_strategy

    def model_dump_json(self, *args, **kwargs):
        self.reassign_ids()
        return super().model_dump_json(*args, **kwargs)
