#!/home/bjorn/anaconda/envs/bjorn/bin/python
# -*- coding: utf-8 -*-
"""
import os
import sys
os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
dna_dir  = u"/home/bjorn/Dropbox/Public/pydna-DNA-assembly//constructs"
dna_dirs = [x[0] for x in os.walk(dna_dir)]
os.environ["pydna_dna_dirs"] = os.pathsep.join(dna_dirs)
os.environ["PYTHONPATH"] = os.environ["pydna_dna_dirs"]
"""

from spyderlib import start_app

start_app.main()
