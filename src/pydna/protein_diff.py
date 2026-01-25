# -*- coding: utf-8 -*-
from Bio.Align import PairwiseAligner
from Bio.Align.substitution_matrices import load

# ---------------------------------------------------------------------
# Input sequences
# ---------------------------------------------------------------------

seq1 = "MKTAYIAKQRQISFVKSHFSRQ"
seq2 = "MKTAYIAKQISFVKSHFSR"

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
        seq = "".join(current["residues"])
        start = current["start"]
        end = current["end"]
        if start == end:
            edits.append(f"Delete {seq} at position {start}")
        else:
            edits.append(f"Delete {seq} at position {start}-{end}")

    elif current["type"] == "ins":
        seq = "".join(current["residues"])
        pos = current["after"]
        edits.append(f"Insert {seq} after position {pos}")

    current = None


# ---------------------------------------------------------------------
# Walk the alignment traceback
# ---------------------------------------------------------------------

p1 = 0  # position in seq1 (0-based)
p2 = 0  # position in seq2 (0-based)

for (s1_start, s1_end), (s2_start, s2_end) in zip(blocks1, blocks2):

    # --- deletions (gap in seq2) ---
    while p1 < s1_start:
        if current and current["type"] == "del" and current["end"] + 1 == p1 + 1:
            current["residues"].append(seq1[p1])
            current["end"] += 1
        else:
            flush()
            current = {
                "type": "del",
                "residues": [seq1[p1]],
                "start": p1 + 1,
                "end": p1 + 1,
            }
        p1 += 1

    # --- insertions (gap in seq1) ---
    while p2 < s2_start:
        if current and current["type"] == "ins" and current["after"] == p1:
            current["residues"].append(seq2[p2])
        else:
            flush()
            current = {
                "type": "ins",
                "residues": [seq2[p2]],
                "after": p1,
            }
        p2 += 1

    # --- aligned region (matches / substitutions) ---
    for i in range(s1_end - s1_start):
        a = seq1[s1_start + i]
        b = seq2[s2_start + i]
        if a != b:
            flush()
            edits.append(f"Substitute {a} → {b} at position {s1_start + i + 1}")

    p1 = s1_end
    p2 = s2_end

# --- tail deletions ---
while p1 < len(seq1):
    if current and current["type"] == "del" and current["end"] + 1 == p1 + 1:
        current["residues"].append(seq1[p1])
        current["end"] += 1
    else:
        flush()
        current = {
            "type": "del",
            "residues": [seq1[p1]],
            "start": p1 + 1,
            "end": p1 + 1,
        }
    p1 += 1

# --- tail insertions ---
while p2 < len(seq2):
    if current and current["type"] == "ins" and current["after"] == p1:
        current["residues"].append(seq2[p2])
    else:
        flush()
        current = {
            "type": "ins",
            "residues": [seq2[p2]],
            "after": p1,
        }
    p2 += 1

flush()

# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

for e in edits:
    print(e)
