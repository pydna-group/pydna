from pathlib import Path
import re
from pydna.dseq import Dseq
from pydna.dseq_old import Dseq as oDseq

text = Path("test_module_dseq.py").read_text(encoding='utf-8')

def repl(m):
    w = m.group(1)
    c = m.group(2)
    t = m.group(3) or ""
    circular = {", circular=True": True, ", circular=False": False, "": False}[t]
    ods = oDseq(w, c, circular=circular)
    nds = Dseq(ods.to_dsiupac(), circular=circular)
    assert ods.to_dsiupac().encode("ascii") == nds._data
    if not repr(ods).strip() == repr(nds):
        pass
        # print(repr(ods))
        # print()
        # print(repr(nds))
        # print("---")
    replacement = f"Dseq(\"{nds._data.decode()}\"{t})"
    hej = eval(replacement)
    assert hej == nds
    return replacement #  m.group(0)

newtext = re.sub("Dseq\((?:\"|')([GATCgatc]+)(?:\"|'), (?:\"|')([GATCgatc]+)(?:\"|')(, (circular=(False|True)))?\)", repl, text)

Path("test_module_dseq_NEW.py").write_text(newtext)
