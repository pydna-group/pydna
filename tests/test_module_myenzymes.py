#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import logging
import os


def test_myenzymes(monkeypatch):
    monkeypatch.setenv("pydna_enzymes", "my_test_enzymes.txt")
    from pydna import myenzymes
    from importlib import reload

    reload(myenzymes)
    assert len(list(myenzymes.myenzymes)) == 621


def test_file_is_folder(monkeypatch, caplog):
    monkeypatch.setenv("pydna_enzymes", "subfolder")
    from pydna import myenzymes
    from importlib import reload

    reload(myenzymes)
    assert len(list(myenzymes.myenzymes)) == 0
    assert len(caplog.records) == 1
    record = next(iter(caplog.records))
    assert record.levelno == logging.WARNING
    if os.name.startswith("nt"):
        assert record.message == "subfolder found, but could not be read."
    else:
        assert record.message == "subfolder is a directory."


def test_file_not_exist(monkeypatch, caplog):
    monkeypatch.setenv("pydna_enzymes", "file_that_does_not_exist.txt")
    from pydna import myenzymes
    from importlib import reload

    reload(myenzymes)
    assert len(list(myenzymes.myenzymes)) == 0
    record = next(iter(caplog.records))
    assert len(caplog.records) == 1
    assert record.levelno == logging.WARNING
    assert record.message == "file_that_does_not_exist.txt not found."


def test_IOError(monkeypatch, caplog):
    from io import StringIO

    monkeypatch.setenv("pydna_enzymes", "my_test_enzymes.txt")

    from unittest.mock import patch

    # from unittest.mock import mock_open

    with patch("pydna.myenzymes.open") as mocked_open:
        mocked_open.side_effect = IOError()
        mocked_open.return_value = StringIO("BamHI")
        from pydna import myenzymes
        from importlib import reload

        reload(myenzymes)

    assert len(list(myenzymes.myenzymes)) == 0
    assert len(caplog.records) == 1
    record = next(iter(caplog.records))
    assert record.levelno == logging.WARNING
    assert record.message == "my_test_enzymes.txt found, but could not be read."


def test_Exception(monkeypatch, caplog):
    from io import StringIO

    monkeypatch.setenv("pydna_enzymes", "my_test_enzymes.txt")

    from unittest.mock import patch

    # from unittest.mock import mock_open

    with patch("pydna.myenzymes.open") as mocked_open:
        mocked_open.side_effect = Exception()
        mocked_open.return_value = StringIO("BamHI")
        from pydna import myenzymes
        from importlib import reload

        reload(myenzymes)

    assert len(list(myenzymes.myenzymes)) == 0
    assert len(caplog.records) == 1
    record = next(iter(caplog.records))
    assert record.levelno == logging.WARNING
    assert "Exception" in record.message


# def test_IOError(monkeypatch, caplog):
#     monkeypatch = pytest.MonkeyPatch()
#     monkeypatch.setenv("pydna_enzymes", "my_test_enzymes.txt")
#     from pydna import myenzymes
#     from importlib import reload
#     reload(myenzymes)
#     from io import StringIO
#     from unittest.mock import patch
#     from unittest.mock import mock_open

#     with patch("myenzymes.open") as mocked_open:
#         mocked_open.side_effect = IOError()
#         mocked_open.return_value = StringIO("BamHI")
#         from pydna import myenzymes
#         from importlib import reload
#         reload(myenzymes)

#     assert len(list(myenzymes.myenzymes)) == 0
#     assert len(caplog.records) == 1
#     record = next(iter(caplog.records))
#     assert record.levelno == logging.WARNING
#     assert record.message == "my_test_enzymes.txt found, but could not be read."

# def test_Exception(monkeypatch, caplog):
#     monkeypatch.setenv("pydna_enzymes", "my_test_enzymes.txt")
#     from pydna import myenzymes
#     from importlib import reload
#     reload(myenzymes)
#     from io import StringIO
#     from pydna import myenzymes
#     from unittest.mock import patch
#     from unittest.mock import mock_open

#     with patch("pydna.myenzymes.open") as mocked_open:
#         mocked_open.side_effect = Exception()
#         mocked_open.return_value = StringIO("BamHI")
#         from pydna import myenzymes
#         from importlib import reload

#         reload(myenzymes)

#     assert len(list(myenzymes.myenzymes)) == 0
#     assert len(caplog.records) == 1
#     record = next(iter(caplog.records))
#     assert record.levelno == logging.WARNING
#     assert "Exception" in record.message


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
