#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
test parse
"""

import pytest


def test_spanning_features():
    pass


#    temp = '''
#
#    LOCUS       MyTemplate                48 bp ds-DNA     circular     04-MAY-2017
#    DEFINITION  .
#    ACCESSION   MyTemplate
#    VERSION     MyTemplate
#    KEYWORDS    .
#    SOURCE
#      ORGANISM  . .
#    COMMENT
#    COMMENT     ApEinfo:methylated:1
#    FEATURES             Location/Qualifiers
#         gene            join(34..48,1..12)
#                         /label=New Feature
#                         /ApEinfo_fwdcolor="cyan"
#                         /ApEinfo_revcolor="green"
#                         /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
#                         width 5 offset 0"
#    ORIGIN
#            1 gctactacac acgtactgac tgcctccaag atagagtcag taaccaca
#    //
#
#    '''
#
#    # feature is
#    #                      gagtcagtaaccacagctactacacac
#
#    tmp = read(temp)
#
#    new_feature = tmp.features[0]
#
#    assert str(new_feature.extract(tmp).seq) == "GAGTCAGTAACCACAGCTACTACACAC"
#
#    new_feature.location.parts[0]
#
#    new_feature.location.parts[1]
#
#    x=tmp._multiply_circular(2)
#
#    assert str(x.features[0].extract(x).seq) == "GAGTCAGTAACCACAGCTACTACACAC"
#
#    temp2 = '''
#
#    LOCUS       MyTemplate                48 bp ds-DNA     circular     04-MAY-2017
#    DEFINITION  .
#    ACCESSION   MyTemplate
#    VERSION     MyTemplate
#    KEYWORDS    .
#    SOURCE
#      ORGANISM  . .
#    COMMENT
#    COMMENT     ApEinfo:methylated:1
#    FEATURES             Location/Qualifiers
#         gene            complement(join(34..48,1..12))
#                         /label=New Feature
#                         /ApEinfo_fwdcolor="cyan"
#                         /ApEinfo_revcolor="green"
#                         /ApEinfo_graphicformat="arrow_data {{0 1 2 0 0 -1} {} 0}
#                         width 5 offset 0"
#    ORIGIN
#            1 gctactacac acgtactgac tgcctccaag atagagtcag taaccaca
#    //
#
#    '''
#
#    tmp2 = read(temp2)
#
#    new_feature2 = tmp2.features[0]
#
#    assert str(new_feature2.extract(tmp2).seq) == "GTGTGTAGTAGCTGTGGTTACTGACTC"
#
#    x2=tmp2._multiply_circular(2)
#
#    assert str(x2.features[0].extract(x2).seq) == "GTGTGTAGTAGCTGTGGTTACTGACTC"

if __name__ == "__main__":
    pytest.cmdline.main([__file__, "-v", "-s"])
