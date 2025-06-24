# -*- coding: utf-8 -*-
import textwrap as _textwrap
import networkx as _nx
from pydna._pretty import pretty_str as _pretty_str
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.utils import rc as _rc


class Contig(_Dseqrecord):
    """This class holds information about a DNA assembly. This class is instantiated by
    the :class:`Assembly` class and is not meant to be used directly.

    """

    def __init__(self, record, *args, graph=None, nodemap=None, **kwargs):
        super().__init__(record, *args, **kwargs)
        self.graph = graph
        self.nodemap = nodemap

    @classmethod
    def from_string(cls, record: str = "", *args, graph=None, nodemap=None, **kwargs):
        obj = super().from_string(record, *args, **kwargs)
        obj.graph = graph
        obj.nodemap = nodemap
        return obj

    @classmethod
    def from_SeqRecord(cls, record, *args, graph=None, nodemap=None, **kwargs):
        obj = super().from_SeqRecord(record, *args, **kwargs)
        obj.graph = graph
        obj.nodemap = nodemap
        return obj

    def __repr__(self):
        return "Contig({}{})".format(
            {True: "-", False: "o"}[not self.circular], len(self)
        )

    def _repr_pretty_(self, p, cycle):
        """returns a short string representation of the object"""
        p.text(
            "Contig({}{})".format({True: "-", False: "o"}[not self.circular], len(self))
        )

    def _repr_html_(self):
        return "<pre>" + self.figure() + "</pre>"

    def reverse_complement(self):
        answer = type(self)(super().reverse_complement())
        g = _nx.DiGraph()
        nm = self.nodemap
        g.add_edges_from(
            [(nm[v], nm[u], d) for u, v, d in list(self.graph.edges(data=True))[::-1]]
        )
        g.add_nodes_from((nm[n], d) for n, d in list(self.graph.nodes(data=True))[::-1])
        for u, v, ed in g.edges(data=True):
            ed["name"] = (
                ed["name"][:-3]
                if ed["name"].endswith("_rc")
                else "{}_rc".format(ed["name"])[:13]
            )
            ed["seq"] = _rc(ed["seq"])
            ln = len(ed["seq"])
            start, stop = ed["piece"].start, ed["piece"].stop
            ed["piece"] = slice(
                ln - stop - g.nodes[u]["length"], ln - start - g.nodes[v]["length"]
            )
            ed["features"] = [f._flip(ln) for f in ed["features"]]
        answer.graph = g
        answer.nodemap = {v: k for k, v in self.nodemap.items()}
        return answer

    rc = reverse_complement

    def detailed_figure(self):
        """Returns a text representation of the assembled fragments.

        Linear:

        ::

            acgatgctatactgCCCCCtgtgctgtgctcta
                               TGTGCTGTGCTCTA
                               tgtgctgtgctctaTTTTTtattctggctgtatc



        Circular:

        ::

            ||||||||||||||
            acgatgctatactgCCCCCtgtgctgtgctcta
                               TGTGCTGTGCTCTA
                               tgtgctgtgctctaTTTTTtattctggctgtatc
                                                  TATTCTGGCTGTATC
                                                  tattctggctgtatcGGGGGtacgatgctatactg
                                                                       ACGATGCTATACTG


        """

        fig = ""
        fragmentposition = 0
        nodeposition = 0
        mylist = []
        for u, v, e in self.graph.edges(data=True):
            nodeposition += e["piece"].stop - e["piece"].start
            fragmentposition -= e["piece"].start
            mylist.append([fragmentposition, e["seq"]])
            mylist.append([nodeposition, v.upper()])
            fragmentposition += e["piece"].stop

        if self.circular:
            edges = list(self.graph.edges(data=True))
            nodeposition = edges[0][2]["piece"].start
            nodelength = len(v)
            mylist = [[nodeposition, "|" * nodelength]] + mylist
        else:
            mylist = mylist[:-1]

        firstpos = -1 * min(0, min(mylist)[0])

        for p, s in mylist:
            fig += "{}{}\n".format(" " * (p + firstpos), s)

        return _pretty_str(fig)

    def figure(self):
        r"""Compact ascii representation of the assembled fragments.

        Each fragment is represented by:

        ::

         Size of common 5' substring|Name and size of DNA fragment|
         Size of common 5' substring

        Linear:

        ::

          frag20| 6
                 \\/
                 /\\
                  6|frag23| 6
                           \\/
                           /\\
                            6|frag14


        Circular:

        ::

          -|2577|61
         |       \\/
         |       /\\
         |       61|5681|98
         |               \\/
         |               /\\
         |               98|2389|557
         |                       \\/
         |                       /\\
         |                       557-
         |                          |
          --------------------------


        """
        nodes = list(self.graph.nodes(data=True))
        edges = list(self.graph.edges(data=True))

        if not self.circular:
            r"""
            frag20| 6
                   \/
                   /\
                    6|frag23| 6
                             \/
                             /\
                              6|frag14
            """

            f = edges[0]

            space2 = len(f[2]["name"])

            fig = ("{name}|{o2:>2}\n" "{space2} \\/\n" "{space2} /\\\n").format(
                name=f[2]["name"], o2=nodes[1][1]["length"], space2=" " * space2
            )
            space = space2  # len(f.name)

            for i, f in enumerate(edges[1:-1]):
                name = "{o1:>2}|{name}|".format(
                    o1=nodes[i + 1][1]["length"], name=f[2]["name"]
                )
                space2 = len(name)

                fig += (
                    "{space} {name}{o2:>2}\n"
                    "{space} {space2}\\/\n"
                    "{space} {space2}/\\\n"
                ).format(
                    name=name,
                    o2=nodes[i + 2][1]["length"],
                    space=" " * space,
                    space2=" " * space2,
                )

                space += space2

            f = edges[-1]
            fig += ("{space} {o1:>2}|{name}").format(
                name=f[2]["name"], o1=nodes[-2][1]["length"], space=" " * (space)
            )

        else:  # circular
            r"""
             -|2577|61
            |       \/
            |       /\
            |       61|5681|98
            |               \/
            |               /\
            |               98|2389|557
            |                       \/
            |                       /\
            |                       557-
            |                          |
             --------------------------
            """

            nodes.append(nodes[0])
            f = edges[0]

            space = len(f[2]["name"]) + 3

            fig = (" -|{name}|{o2:>2}\n" "|{space}\\/\n" "|{space}/\\\n").format(
                name=f[2]["name"], o2=nodes[1][1]["length"], space=" " * space
            )

            for i, f in enumerate(edges[1:]):
                name = "{o1:>2}|{name}|".format(
                    o1=nodes[i + 1][1]["length"], name=f[2]["name"]
                )
                space2 = len(name)
                fig += (
                    "|{space}{name}{o2:>2}\n"
                    "|{space}{space2}\\/\n"
                    "|{space}{space2}/\\\n"
                ).format(
                    o2=nodes[i + 2][1]["length"],
                    name=name,
                    space=" " * space,
                    space2=" " * space2,
                )
                space += space2

            fig += "|{space}{o1:>2}-\n".format(
                space=" " * (space), o1=nodes[0][1]["length"]
            )
            fig += "|{space}   |\n".format(space=" " * (space))
            fig += " {space}".format(space="-" * (space + 3))
        return _pretty_str(_textwrap.dedent(fig))
