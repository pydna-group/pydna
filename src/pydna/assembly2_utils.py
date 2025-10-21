# -*- coding: utf-8 -*-
"""
Functions that are used in the assembly2 module but not part of the class, and that may be
useful in other parts of the package. Separated to avoid circular imports.
"""

from pydna.types import EdgeRepresentationAssembly, SubFragmentRepresentationAssembly
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.utils import location_boundaries as _location_boundaries


def is_sublist(sublist: list, my_list: list, my_list_is_cyclic: bool = False) -> bool:
    """Returns True if argument sublist is a sublist of argument my_list (can be treated as cyclic), False otherwise.

    Examples
    --------
    >>> is_sublist([1, 2], [1, 2, 3], False)
    True
    >>> is_sublist([1, 2], [1, 3, 2], False)
    False

    # See the case here for cyclic lists
    >>> is_sublist([3, 1], [1, 2, 3], False)
    False
    >>> is_sublist([3, 1], [1, 2, 3], True)
    True
    """
    n = len(sublist)
    if my_list_is_cyclic:
        my_list = my_list + my_list
    for i in range(len(my_list) - n + 1):
        # Just in case tuples were passed
        if list(my_list[i : i + n]) == list(sublist):
            return True
    return False


def reverse_complement_assembly(
    assembly: EdgeRepresentationAssembly, fragments: list[_Dseqrecord]
) -> EdgeRepresentationAssembly:
    """Complement an assembly, i.e. reverse the order of the fragments and the orientation of the overlaps."""
    new_assembly = list()
    for u, v, locu, locv in assembly:
        f_u = fragments[abs(u) - 1]
        f_v = fragments[abs(v) - 1]
        new_assembly.append((-v, -u, locv._flip(len(f_v)), locu._flip(len(f_u))))
    return new_assembly[::-1]


def filter_linear_subassemblies(
    linear_assemblies: list[EdgeRepresentationAssembly],
    circular_assemblies: list[EdgeRepresentationAssembly],
    fragments: list[_Dseqrecord],
) -> list[EdgeRepresentationAssembly]:
    """Remove linear assemblies which are sub-assemblies of circular assemblies"""
    all_circular_assemblies = circular_assemblies + [
        reverse_complement_assembly(c, fragments) for c in circular_assemblies
    ]
    filtered_assemblies = [
        assem
        for assem in linear_assemblies
        if not any(is_sublist(assem, c, True) for c in all_circular_assemblies)
    ]
    # I don't think the line below is necessary, but just in case
    # filtered_assemblies = [l for l in filtered_assemblies if not any(is_sublist(reverse_complement_assembly(l, fragments), c, True) for c in all_circular_assemblies)]
    return filtered_assemblies


def remove_subassemblies(
    assemblies: list[EdgeRepresentationAssembly],
) -> list[EdgeRepresentationAssembly]:
    """Filter out subassemblies, i.e. assemblies that are contained within another assembly.

    For example:
        [(1, 2, '1[8:14]:2[1:7]'), (2, 3, '2[10:17]:3[1:8]')]
        [(1, 2, '1[8:14]:2[1:7]')]
    The second one is a subassembly of the first one.
    """

    # Sort by length, longest first
    assemblies = sorted(assemblies, key=len, reverse=True)

    filtered_assemblies = list()
    for assembly in assemblies:
        # Check if this assembly is a subassembly of any of the assemblies we have already found
        if not any(is_sublist(assembly, a) for a in filtered_assemblies):
            filtered_assemblies.append(assembly)

    return filtered_assemblies


def assembly2str(assembly: EdgeRepresentationAssembly) -> str:
    """Convert an assembly to a string representation, for example:
    ((1, 2, [8:14], [1:7]),(2, 3, [10:17], [1:8]))
    becomes:
    ('1[8:14]:2[1:7]', '2[10:17]:3[1:8]')

    The reason for this is that by default, a feature '[8:14]' when present in a tuple
    is printed to the console as ``SimpleLocation(ExactPosition(8), ExactPosition(14), strand=1)`` (very long).
    """
    return str(tuple(f"{u}{lu}:{v}{lv}" for u, v, lu, lv in assembly))


def assembly2str_tuple(assembly: EdgeRepresentationAssembly) -> str:
    """Convert an assembly to a string representation, like
    ((1, 2, [8:14], [1:7]),(2, 3, [10:17], [1:8]))
    """
    return str(tuple((u, v, str(lu), str(lv)) for u, v, lu, lv in assembly))


def assembly_has_mismatches(
    fragments: list[_Dseqrecord], assembly: EdgeRepresentationAssembly
) -> bool:
    """Check if an assembly has mismatches. This should never happen and if so it returns an error."""
    for u, v, loc_u, loc_v in assembly:
        seq_u = fragments[u - 1] if u > 0 else fragments[-u - 1].reverse_complement()
        seq_v = fragments[v - 1] if v > 0 else fragments[-v - 1].reverse_complement()
        # TODO: Check issue where extraction failed, and whether it would give problems here
        if (
            str(loc_u.extract(seq_u).seq).upper()
            != str(loc_v.extract(seq_v).seq).upper()
        ):
            return True
    return False


def assembly_is_circular(
    assembly: EdgeRepresentationAssembly, fragments: list[_Dseqrecord]
) -> bool:
    """
    Based on the topology of the locations of an assembly, determine if it is circular.
    This does not work for insertion assemblies, that's why assemble takes the optional argument is_insertion.
    """
    if assembly[0][0] != assembly[-1][1]:
        return False
    elif (
        isinstance(fragments[abs(assembly[0][0]) - 1], _Dseqrecord)
        and fragments[abs(assembly[0][0]) - 1].circular
    ):
        return True
    else:
        return (
            _location_boundaries(assembly[0][2])[0]
            > _location_boundaries(assembly[-1][3])[0]
        )


def edge_representation2subfragment_representation(
    assembly: EdgeRepresentationAssembly, is_circular: bool
) -> SubFragmentRepresentationAssembly:
    """
    Turn this kind of edge representation fragment 1, fragment 2, right edge on 1, left edge on 2
    a = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b', 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    Into this: fragment 1, left edge on 1, right edge on 1
    b = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    """

    if is_circular:
        temp = list(assembly[-1:]) + list(assembly)
    else:
        temp = (
            [(None, assembly[0][0], None, None)]
            + list(assembly)
            + [(assembly[-1][1], None, None, None)]
        )
    edge_pairs = zip(temp, temp[1:])
    subfragment_representation = list()
    for (_u1, v1, _, start_location), (_u2, _v2, end_location, _) in edge_pairs:
        subfragment_representation.append((v1, start_location, end_location))

    return tuple(subfragment_representation)


def subfragment_representation2edge_representation(
    assembly: SubFragmentRepresentationAssembly, is_circular: bool
) -> EdgeRepresentationAssembly:
    """
    Turn this kind of subfragment representation fragment 1, left edge on 1, right edge on 1
    a = [(1, 'loc1c', 'loc1a'), (2, 'loc2a', 'loc2b'), (3, 'loc3b', 'loc3c')]
    Into this: fragment 1, fragment 2, right edge on 1, left edge on 2
    b = [(1, 2, 'loc1a', 'loc2a'), (2, 3, 'loc2b' 'loc3b'), (3, 1, 'loc3c', 'loc1c')]
    """

    edge_representation = []

    # Iterate through the assembly pairwise to create the edge representation
    for i in range(len(assembly) - 1):
        frag1, left1, right1 = assembly[i]
        frag2, left2, right2 = assembly[i + 1]
        # Create the edge between the current and next fragment
        edge_representation.append((frag1, frag2, right1, left2))

    if is_circular:
        # Add the edge from the last fragment back to the first
        frag_last, left_last, right_last = assembly[-1]
        frag_first, left_first, right_first = assembly[0]
        edge_representation.append((frag_last, frag_first, right_last, left_first))

    return tuple(edge_representation)
