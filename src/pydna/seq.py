#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A subclass of the Biopython SeqRecord class.

Has a number of extra methods and uses
the :class:`pydna._pretty_str.pretty_str` class instread of str for a
nicer output in the IPython shell.
"""

# from pydna.codon import weights as _weights
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pydna.codon import rare_codons as _rare_codons
from pydna.codon import start as _start
from pydna.codon import stop as _stop
from pydna.codon import n_end as _n_end
from seguid import lsseguid as _lsseguid
from pydna.utils import rc as _rc

from Bio.SeqUtils import seq3 as _seq3
from Bio.SeqUtils import gc_fraction as _GC
import re as _re
from Bio.Seq import Seq as _Seq
from pydna._pretty import PrettyTable as _PrettyTable

from typing import List as _List, Optional as _Optional, Tuple as _Tuple

# import logging as _logging

# _module_logger = _logging.getLogger("pydna." + __name__)


class Seq(_Seq):
    """docstring."""

    def translate(
        self,
        *args,
        stop_symbol: str = "*",
        to_stop: bool = False,
        cds: bool = False,
        gap: str = "-",
        **kwargs,
    ) -> "ProteinSeq":
        """Translate.."""
        p = super().translate(
            *args, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap, **kwargs
        )
        return ProteinSeq(p._data)

    def gc(self) -> float:
        """Return GC content."""
        return round(_GC(self._data.upper().decode("ASCII")), 3)

    def cai(self, organism: str = "sce") -> float:
        """docstring."""
        from pydna.utils import cai as _cai

        return _cai(self._data.upper().decode("ASCII"), organism=organism)

    def rarecodons(self, organism: str = "sce") -> _List[slice]:
        """docstring."""
        rare = _rare_codons[organism]
        s = self._data.upper().decode("ASCII")
        slices: _List[slice] = []
        for i in range(0, len(self) // 3):
            x, y = i * 3, i * 3 + 3
            trip = s[x:y]
            if trip in rare:
                slices.append(slice(x, y, 1))
        return slices

    def startcodon(self, organism: str = "sce") -> _Optional[float]:
        """docstring."""
        return _start[organism].get(self._data.upper().decode("ASCII")[:3])

    def stopcodon(self, organism: str = "sce") -> _Optional[float]:
        """docstring."""
        return _stop[organism].get(self._data.upper().decode("ASCII")[-3:])

    def express(self, organism: str = "sce") -> _PrettyTable:
        """docstring."""
        x = _PrettyTable(
            ["cds", "len", "cai", "gc", "sta", "stp", "n-end"]
            + _rare_codons[organism]
            + ["rare"]
        )
        val = []

        val.append(
            f"{self._data.upper().decode('ASCII')[:3]}..."
            f"{self._data.upper().decode('ASCII')[-3:]}"
        )
        val.append(len(self) / 3)
        val.append(self.cai(organism))
        val.append(self.gc())
        val.append(self.startcodon())
        val.append(self.stopcodon())
        val.append(
            _n_end[organism].get(_seq3(self[3:6].translate())),
        )
        s = self._data.upper().decode("ASCII")
        trps = [s[i * 3 : i * 3 + 3] for i in range(0, len(s) // 3)]
        tot = 0
        for cdn in _rare_codons[organism]:
            cnt = trps.count(cdn)
            tot += cnt
            val.append(cnt)
        val.append(round(tot / len(trps), 3))
        x.add_row(val)
        return x

    def orfs2(self, minsize: int = 30) -> _List[str]:
        """docstring."""
        orf = _re.compile(
            f"ATG(?:...){{{minsize},}}?(?:TAG|TAA|TGA)", flags=_re.IGNORECASE
        )
        start = 0
        matches: _List[slice] = []
        s = self._data.decode("ASCII")

        while True:
            match = orf.search(s, pos=start)
            if match:
                matches.append(slice(match.start(), match.end()))
                start = match.start() + 1
            else:
                break
        return sorted([self[sl] for sl in matches], key=len, reverse=True)

    def orfs(self, minsize: int = 100) -> _List[_Tuple[int, int]]:
        dna = self._data.decode("ASCII")
        from pydna.utils import three_frame_orfs

        return [(x, y) for frame, x, y in three_frame_orfs(dna, limit=minsize)]

    def seguid(self) -> str:
        """Url safe SEGUID [#]_ for the sequence.

        This checksum is the same as seguid but with base64.urlsafe
        encoding instead of the normal base64. This means that
        the characters + and / are replaced with - and _ so that
        the checksum can be part of a URL.

        Examples
        --------
        >>> from pydna.seq import Seq
        >>> a = Seq("aa")
        >>> a.seguid()
        'lsseguid=gBw0Jp907Tg_yX3jNgS4qQWttjU'

        References
        ----------
        .. [#] http://wiki.christophchamp.com/index.php/SEGUID
        """
        return _lsseguid(self._data.decode("utf8").upper(), alphabet="{DNA-extended}")

    def __getitem__(self, key):
        result = super().__getitem__(key)
        try:
            result.__class__ = self.__class__
        except TypeError:
            pass
        return result

    def reverse_complement(self):
        return self.__class__(_rc(self._data))

    rc = reverse_complement


class ProteinSeq(_Seq):
    """docstring."""

    def translate(self):
        raise NotImplementedError("Not defined for protein.")

    def complement(self):
        raise NotImplementedError("Not defined for protein.")

    def complement_rna(self):
        raise NotImplementedError("Not defined for protein.")

    def reverse_complement(self):
        raise NotImplementedError("Not defined for protein.")

    rc = reverse_complement

    def reverse_complement_rna(self):
        raise NotImplementedError("Not defined for protein.")

    def transcribe(self):
        raise NotImplementedError("Not defined for protein.")

    def back_transcribe(self):
        raise NotImplementedError("Not defined for protein.")

    def seguid(self) -> str:
        """Url safe SEGUID [#]_ for the sequence.

        This checksum is the same as seguid but with base64.urlsafe
        encoding instead of the normal base64. This means that
        the characters + and / are replaced with - and _ so that
        the checksum can be part of a URL.

        Examples
        --------
        >>> from pydna.seq import ProteinSeq
        >>> a = ProteinSeq("aa")
        >>> a.seguid()
        'lsseguid=gBw0Jp907Tg_yX3jNgS4qQWttjU'

        References
        ----------
        .. [#] http://wiki.christophchamp.com/index.php/SEGUID
        """
        return _lsseguid(
            self._data.decode("utf8").upper(), alphabet="{protein-extended}"
        )

    def __getitem__(self, key):
        result = super().__getitem__(key)
        try:
            result.__class__ = self.__class__
        except TypeError:
            pass
        return result

    def _pa(self) -> ProteinAnalysis:
        # breakpoint()
        return ProteinAnalysis(self._data.decode("ascii"))

    def molecular_weight(self) -> float:
        return self._pa().molecular_weight()

    def pI(self) -> float:
        return self._pa().isoelectric_point()

    def instability_index(self) -> float:
        """
        Instability index according to Guruprasad et al.

        Value above 40 means the protein is has a short half life.

        Guruprasad K., Reddy B.V.B., Pandit M.W. Protein Engineering 4:155-161(1990).
        """
        return self._pa().instability_index()
