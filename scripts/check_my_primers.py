#!/home/bjorn/anaconda3/envs/bjorn36/bin/python
# -*- coding: utf-8 -*-


from pydna.myprimers import list_primers

unique = set(str(p.seq).lower() for p in list_primers)

print("number of primers", len(list_primers))
print("unique primers", len(unique))
print("no seq (n)", len([p for p in list_primers if set(p.seq.lower()) == set("n")]))

for i, p in enumerate(list_primers):
    if not p.name.startswith(str(i)):
        print(i, p.format("tab"))

print("names checked for primer number for ", i + 1, "primers")

from collections import defaultdict

dct = defaultdict(list)

for u in unique:
    for p in list_primers:
        if set(u) == set("n") or set(p.seq.lower()) == set("n"):
            continue
        if u == str(p.seq).lower():
            dct[u].append(p.name)


for seq, names in dct.items():
    if len(names) > 1:
        print("\n".join(names))
        print(seq)
        print()
