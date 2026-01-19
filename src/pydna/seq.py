#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A subclass of Biopython Bio.Seq.Seq

Has a number of extra methods and uses
the :class:`pydna._pretty_str.pretty_str` class instread of str for a
nicer output in the IPython shell.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pydna.codon import rare_codons
from pydna.codon import start as _start
from pydna.codon import stop as _stop
from pydna.codon import n_end as _n_end
from seguid import lsseguid
from pydna.utils import rc

from Bio.SeqUtils import seq3
from Bio.SeqUtils import gc_fraction
import re
from Bio.Seq import Seq as _Seq
from pydna._pretty import PrettyTable

from typing import List, Optional, Tuple


class Seq(_Seq):
    """docstring."""

    # @property
    # def full_sequence(self):
    #     return self

    # def translate(
    #     self,
    #     *args,
    #     stop_symbol: str = "*",
    #     to_stop: bool = False,
    #     cds: bool = False,
    #     gap: str = "-",
    #     **kwargs,
    # ) -> "ProteinSeq":
    #     """Translate.."""
    #     p = super().translate(
    #         *args, stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap, **kwargs
    #     )
    #     return ProteinSeq(p._data)

    def translate(
        self,
        table: [str, int] = "Standard",
        stop_symbol: [str] = "*",
        to_stop: bool = False,
        cds: bool = False,
        gap: str = "-",
    ) -> _Seq:

        # TODO: is this method needed?
        """
        Translate into protein.

        The table argument is the name of a codon table (string). These names
        can be for example "Standard" or "Alternative Yeast Nuclear" for the
        yeast CUG clade where the CUG codon is translated as serine instead
        of the standard leucine.

        Over forty translation tables are available from the BioPython
        Bio.Data.CodonTable module. Look at the keys of the dictionary
        ´CodonTable.ambiguous_generic_by_name´.
        These are based on tables in this file provided by NCBI:

        https://ftp.ncbi.nlm.nih.gov/entrez/misc/data/gc.prt

        Standard table

          |  T      |  C      |  A      |  G      |
        --+---------+---------+---------+---------+--
        T | TTT F   | TCT S   | TAT Y   | TGT C   | T
        T | TTC F   | TCC S   | TAC Y   | TGC C   | C
        T | TTA L   | TCA S   | TAA Stop| TGA Stop| A
        T | TTG L(s)| TCG S   | TAG Stop| TGG W   | G
        --+---------+---------+---------+---------+--
        C | CTT L   | CCT P   | CAT H   | CGT R   | T
        C | CTC L   | CCC P   | CAC H   | CGC R   | C
        C | CTA L   | CCA P   | CAA Q   | CGA R   | A
        C | CTG L(s)| CCG P   | CAG Q   | CGG R   | G
        --+---------+---------+---------+---------+--
        A | ATT I   | ACT T   | AAT N   | AGT S   | T
        A | ATC I   | ACC T   | AAC N   | AGC S   | C
        A | ATA I   | ACA T   | AAA K   | AGA R   | A
        A | ATG M(s)| ACG T   | AAG K   | AGG R   | G
        --+---------+---------+---------+---------+--
        G | GTT V   | GCT A   | GAT D   | GGT G   | T
        G | GTC V   | GCC A   | GAC D   | GGC G   | C
        G | GTA V   | GCA A   | GAA E   | GGA G   | A
        G | GTG V   | GCG A   | GAG E   | GGG G   | G
        --+---------+---------+---------+---------+--


        Parameters
        ----------
        table : [str, int], optional
            The default is "Standard". Can be a table id integer, see here for table
            numbering https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
        stop_symbol : [str], optional
            The default is "*". Single character string to indicate translation stop.
        to_stop : bool, optional
            The default is False. True means that translation terminates at the first
              in frame stop codon. False translates to the end.
        cds : bool, optional
            The default is False. If True, checks that the sequence starts with a
            valid alternative start codon sequence length is a multiple of three, and
            that there is a single in frame stop codon at the end. If these tests fail,
            an exception is raised.
        gap : str, optional
            The default is "-".

        Returns
        -------
        Bio.Seq.Seq
            A Biopython Seq object with the translated amino acid code.

        """

        p = _Seq(self._data).translate(
            stop_symbol=stop_symbol, to_stop=to_stop, cds=cds, gap=gap
        )
        return ProteinSeq(p._data)

    def transcribe(self) -> _Seq:
        """
        Transcribe a DNA sequence into RNA and return the RNA sequence
        as a new Seq object.

        """
        return Seq(_Seq(self._data).transcribe()._data)

    def gc(self) -> float:
        """Return GC content."""
        return round(gc_fraction(self._data.upper().decode("ASCII")), 3)

    def cai(self, organism: str = "sce") -> float:
        """docstring."""
        from pydna.utils import cai as _cai

        return _cai(self._data.upper().decode("ASCII"), organism=organism)

    def rarecodons(self, organism: str = "sce") -> List[slice]:
        """docstring."""
        rare = rare_codons[organism]
        s = self._data.upper().decode("ASCII")
        slices: List[slice] = []
        for i in range(0, len(self) // 3):
            x, y = i * 3, i * 3 + 3
            trip = s[x:y]
            if trip in rare:
                slices.append(slice(x, y, 1))
        return slices

    def startcodon(self, organism: str = "sce") -> Optional[float]:
        """docstring."""
        return _start[organism].get(self._data.upper().decode("ASCII")[:3])

    def stopcodon(self, organism: str = "sce") -> Optional[float]:
        """docstring."""
        return _stop[organism].get(self._data.upper().decode("ASCII")[-3:])

    def express(self, organism: str = "sce") -> PrettyTable:
        """docstring."""
        x = PrettyTable(
            ["cds", "len", "cai", "gc", "sta", "stp", "n-end"]
            + rare_codons[organism]
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
            _n_end[organism].get(seq3(self[3:6].translate())),
        )
        s = self._data.upper().decode("ASCII")
        trps = [s[i * 3 : i * 3 + 3] for i in range(0, len(s) // 3)]
        tot = 0
        for cdn in rare_codons[organism]:
            cnt = trps.count(cdn)
            tot += cnt
            val.append(cnt)
        val.append(round(tot / len(trps), 3))
        x.add_row(val)
        return x

    def orfs2(self, minsize: int = 30) -> List[str]:
        """docstring."""
        orf = re.compile(
            f"ATG(?:...){{{minsize},}}?(?:TAG|TAA|TGA)", flags=re.IGNORECASE
        )
        start = 0
        matches: List[slice] = []
        s = self._data.decode("ASCII")

        while True:
            match = orf.search(s, pos=start)
            if match:
                matches.append(slice(match.start(), match.end()))
                start = match.start() + 1
            else:
                break
        return sorted([self[sl] for sl in matches], key=len, reverse=True)

    def orfs(self, minsize: int = 100) -> List[Tuple[int, int]]:
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
        return lsseguid(
            self._data.decode("ascii").upper(), alphabet="{DNA-extended},AU"
        )

    # def __getitem__(self, key):
    #     result = super().__getitem__(key)
    #     try:
    #         result.__class__ = self.__class__
    #     except TypeError:
    #         pass
    #     return result

    def reverse_complement(self):
        return self.__class__(rc(self._data))

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
        return lsseguid(
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
