# -*- coding: utf-8 -*-
import textwrap
import networkx as nx
from pydna._pretty import pretty_str as ps
from pydna.dseqrecord import Dseqrecord
from pydna.utils import rc
import numpy as np


class Contig(Dseqrecord):
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
        g = nx.DiGraph()
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
            ed["seq"] = rc(ed["seq"])
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

        return ps(fig)

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
        return ps(textwrap.dedent(fig))

    def figure_mpl(self):
        """
        Graphic representation of the assembly.

        Returns
        -------
        matplotlib.figure.Figure
            A representation of a linear or culrcular assembly.

        """
        # lazy imports in case matplotlib is not installed
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches

        plt.ioff()  # Disable interactive mode, otherwise two plots are shown in Spyder.
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.isinteractive.html#matplotlib.pyplot.isinteractive

        def pick_n_colors(n, cmap_name="tab20"):
            cmap = plt.get_cmap(cmap_name)
            return [cmap(i / n) for i in range(n)]

        fig, ax = plt.subplots()
        edges = list(self.graph.edges(data=True))
        colors = pick_n_colors(len(edges))

        if self.circular:
            # Circle parameters for Circular assembly
            center = 0, 0
            outer_radius = 1.5  # fragments on the outer lane
            middle_radius = 1.3  # fragments on the inner lane
            small_radius = 1.1  # odd number of fragments require an extra radius
            arc_width = 0.1  # Arc thickness

            circle = len(self)  # The circle has the length of the assembly
            radii = [outer_radius, middle_radius] * (
                len(edges) // 2
            )  # radii alternates, starting with outer.

            if len(edges) % 2 != 0:  # last fragment get a smaller radius
                radii.append(small_radius)

            assert (
                len(colors) == len(radii) == len(edges)
            )  # One color and one radius for each edge.

            # The recombination between last and first fragments
            # end at the origin (twelve o'clock).
            start = 0 - len(edges[0][0])

            for edge, radius, color in zip(edges, radii, colors):

                node1, node2, meta = edge
                slc = meta["piece"]
                extra = len(node2)
                # slc contain the first but not the second node, so add extra to the length
                length = slc.stop - slc.start + extra

                theta1 = 90.0 - 360.0 / circle * start
                theta2 = 90.0 - 360.0 / circle * (start + length)

                # Create arc
                arc_patch = mpatches.Wedge(
                    center=center,
                    r=radius,
                    theta1=theta2,
                    theta2=theta1,
                    width=arc_width,
                    edgecolor=color,
                    facecolor=(1, 1, 1, 0),
                    linewidth=1,
                )
                ax.add_patch(arc_patch)

                # Compute label position slightly outside the arc
                mid_angle = (theta1 + theta2) / 2
                rad = np.deg2rad(mid_angle)
                label_radius = radius + arc_width + 0.1  # place outside the arc
                x = label_radius * np.cos(rad)
                y = label_radius * np.sin(rad)

                # Choose alignment based on angle
                ha = "left" if np.cos(rad) >= 0 else "right"
                va = "center"

                ax.text(x, y, meta["name"], ha=ha, va=va, fontsize=10)

                start += length - len(node2)
            ax.axis("off")
            ax.set_aspect("equal")
            ax.set_xlim(-1.6, 1.6)  # This should be enough, but not extensively tested.
            ax.set_ylim(-1.6, 1.6)

        else:  # Linear assembly
            import itertools  # 3131 bp

            unit = len(self) / 50
            upper = 4 * unit
            lower = 1 * unit
            height = 1 * unit
            x = 0

            for edge, y, color in zip(edges, itertools.cycle((lower, upper)), colors):
                node1, node2, metadict = edge
                slc = metadict["piece"]
                # slc contain the first but not the second node, so add extra to the length if not begin or end.
                extra = len(node2) if node2 not in ("begin", "end") else 0
                length = slc.stop - slc.start + extra
                box = mpatches.FancyBboxPatch(
                    (x, y),
                    length,
                    height,
                    linewidth=1,
                    boxstyle="round",
                    edgecolor=color,
                    facecolor=(1, 1, 1, 0),
                )
                ax.add_patch(box)
                ax.text(
                    x + length / 2,
                    y + height * 2 if y == upper else y - height * 2,
                    metadict["name"],
                    ha="center",
                    va="center",
                    fontsize=10,
                )
                x += length - len(node2)
            ax.axis("off")
            ax.set_aspect("equal")
            ax.set_xlim(-1, len(self) + 1)
            ax.set_ylim(-height, height * 2 + upper)
        return fig

    # FIXME: This code uses plotly, but I see no reason for it at this point.
    # def figure_plotly(self):
    #     import plotly.graph_objects as go
    #     import numpy as np

    #     circ = len(self)
    #     arcs = list(self.graph.edges(data=True))

    #     # Radii setup
    #     small_radius = 1.1
    #     middle_radius = 1.3
    #     outer_radius = 1.5
    #     arc_width = 0.1

    #     radii = [outer_radius, middle_radius] * (len(arcs) // 2)
    #     if len(arcs) % 2 != 0:
    #         radii.append(small_radius)

    #     fig = go.Figure()
    #     start = 0 - len(arcs[0][0])

    #     for (node1, node2, meta), radius in zip(arcs, radii):
    #         slc = meta["piece"]
    #         length = slc.stop - slc.start + len(node1)

    #         theta1 = 90.0 - 360.0 / circ * start
    #         theta2 = 90.0 - 360.0 / circ * (start + length)

    #         # Generate arc points
    #         theta = np.linspace(theta1, theta2, 50)
    #         theta_rev = theta[::-1]

    #         r_outer = np.full_like(theta, radius)
    #         r_inner = np.full_like(theta_rev, radius - arc_width)

    #         r = np.concatenate([r_outer, r_inner])
    #         t = np.concatenate([theta, theta_rev])

    #         fig.add_trace(
    #             go.Scatterpolar(
    #                 r=r,
    #                 theta=t,
    #                 fill="toself",
    #                 mode="lines",
    #                 line_color="rgba(0,100,200,0.6)",
    #                 hoverinfo="text",
    #                 text=meta["name"],
    #                 name=meta["name"],
    #             )
    #         )

    #         start += length - len(node2)

    #     fig.update_layout(
    #         polar=dict(radialaxis=dict(visible=False), angularaxis=dict(visible=False)),
    #         showlegend=False,
    #     )

    #     fig.show("browser")
