# -*- coding: utf-8 -*-
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load

# ---------------------------------------------------------------------
# Input sequences
# ---------------------------------------------------------------------

seq1 = "MKTAYIAKKKKKISFVKSHFSR"
seq2 = "MKTAYIAKQRQISFVKSHFSRQ"
seq1 = "MKTAYIAKQRQISFVKSHFSRQ"
seq2 = "MKTAYIAKQISFVKSHFSR"
seq2 = "MKTAYIAKQRQISFVKSHFSRQ"
seq1 = "MKTAYIAKQISFVKSHFSR"
# ---------------------------------------------------------------------
# Alignment setup
# ---------------------------------------------------------------------

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.substitution_matrix = load("BLOSUM62")
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5

alignment = aligner.align(seq1, seq2)[0]

print(f"Score: {alignment.score}")
print(alignment)

blocks1, blocks2 = alignment.aligned

# ---------------------------------------------------------------------
# Edit accumulation helpers
# ---------------------------------------------------------------------

edits = []
current = None  # open edit block


def flush():
    """Emit the currently accumulated edit (if any)."""
    global current
    if current is None:
        return

    if current["type"] == "del":
        seq = "".join(current["from"])
        s, e = current["start"], current["end"]
        if s == e:
            edits.append(f"Delete {seq} at position {s}")
        else:
            edits.append(f"Delete {seq} at position {s}-{e}")

    elif current["type"] == "ins":
        seq = "".join(current["to"])
        edits.append(f"Insert {seq} after position {current['after']}")

    elif current["type"] == "sub":
        frm = "".join(current["from"])
        to = "".join(current["to"])
        s, e = current["start"], current["end"]
        if s == e:
            edits.append(f"Substitute {frm} → {to} at position {s}")
        else:
            edits.append(f"Substitute {frm} → {to} from position {s} to {e}")

    current = None


# ---------------------------------------------------------------------
# Walk the alignment traceback
# ---------------------------------------------------------------------

p1 = 0  # position in seq1 (0-based)
p2 = 0  # position in seq2 (0-based)

for (s1_start, s1_end), (s2_start, s2_end) in zip(blocks1, blocks2):

    # --- deletions ---
    while p1 < s1_start:
        if current and current["type"] == "del" and current["end"] + 1 == p1 + 1:
            current["from"].append(seq1[p1])
            current["end"] += 1
        else:
            flush()
            current = {
                "type": "del",
                "from": [seq1[p1]],
                "start": p1 + 1,
                "end": p1 + 1,
            }
        p1 += 1

    # --- insertions ---
    while p2 < s2_start:
        if current and current["type"] == "ins" and current["after"] == p1:
            current["to"].append(seq2[p2])
        else:
            flush()
            current = {
                "type": "ins",
                "to": [seq2[p2]],
                "after": p1,
            }
        p2 += 1

    # --- aligned region (matches / substitutions) ---
    for i in range(s1_end - s1_start):
        a = seq1[s1_start + i]
        b = seq2[s2_start + i]
        pos = s1_start + i + 1

        if a == b:
            flush()
            continue

        if current and current["type"] == "sub" and current["end"] + 1 == pos:
            current["from"].append(a)
            current["to"].append(b)
            current["end"] += 1
        else:
            flush()
            current = {
                "type": "sub",
                "from": [a],
                "to": [b],
                "start": pos,
                "end": pos,
            }

    p1 = s1_end
    p2 = s2_end

# --- tail deletions ---
while p1 < len(seq1):
    if current and current["type"] == "del" and current["end"] + 1 == p1 + 1:
        current["from"].append(seq1[p1])
        current["end"] += 1
    else:
        flush()
        current = {
            "type": "del",
            "from": [seq1[p1]],
            "start": p1 + 1,
            "end": p1 + 1,
        }
    p1 += 1

# --- tail insertions ---
while p2 < len(seq2):
    if current and current["type"] == "ins" and current["after"] == p1:
        current["to"].append(seq2[p2])
    else:
        flush()
        current = {
            "type": "ins",
            "to": [seq2[p2]],
            "after": p1,
        }
    p2 += 1

flush()

# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

for e in edits:
    print(e)
