import importlib

import pytest


shims = [
    ("dseq", "Dseq"),
    ("seq", "Seq"),
    ("dseqrecord", "Dseqrecord"),
    ("seqrecord", "SeqRecord"),
    ("alphabet", "get_parts"),
    ("types", "DseqType"),
]


@pytest.mark.parametrize("name,attr", shims)
def test_legacy_shim_forwards_to_core(name, attr):

    legacy = importlib.import_module(f"pydna.legacy.{name}")
    core = importlib.import_module(f"pydna.core.{name}")

    assert getattr(legacy, attr) is getattr(core, attr)


@pytest.mark.parametrize("name,attr", shims)
def test_legacy_shim_dir(name, attr):

    legacy = importlib.import_module(f"pydna.legacy.{name}")

    assert attr in dir(legacy)


@pytest.mark.parametrize("name,attr", shims)
def test_legacy_shim_missing_attribute(name, attr):

    legacy = importlib.import_module(f"pydna.legacy.{name}")

    with pytest.raises(AttributeError):
        legacy.this_attribute_does_not_exist
