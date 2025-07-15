#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import tempfile
import os
import shutil

from unittest import mock

import textwrap


# def test_makedirs_fail(monkeypatch, caplog):

#     from unittest.mock import patch

#     with patch("pydna._os.path.isdir") as pid, patch("pydna._os.makedirs") as pmd:
#         pmd.side_effect = OSError()
#         pid.return_value = True
#         with pytest.raises(OSError):
#             import pydna
#             from importlib import reload

#             reload(pydna)
#         # assert pmd.called == True

#     with patch("pydna._os.path.isdir") as pid, patch("pydna._os.makedirs") as pmd:
#         pmd.side_effect = IOError()  # OSError()
#         pid.return_value = False
#         with pytest.raises(IOError):
#             import pydna
#             from importlib import reload

#             reload(pydna)


# def test_default_env(monkeypatch):
#     pydna_base_dir = os.path.join(tempfile.gettempdir(), "pydna_test")

#     try:
#         shutil.rmtree(pydna_base_dir)
#     except FileNotFoundError:
#         pass

#     def pydna_config_dir_name(x):
#         return pydna_base_dir

#     def pydna_data_dir_name(x):
#         return pydna_base_dir

#     def pydna_log_dir_name(x):
#         return pydna_base_dir

#     monkeypatch.delenv("pydna_loglevel", raising=False)
#     monkeypatch.delenv("pydna_email", raising=False)
#     monkeypatch.delenv("pydna_cached_funcs", raising=False)
#     monkeypatch.delenv("pydna_ape", raising=False)
#     monkeypatch.delenv("pydna_primers", raising=False)
#     monkeypatch.delenv("pydna_enzymes", raising=False)
#     monkeypatch.delenv("pydna_config_dir", raising=False)
#     monkeypatch.delenv("pydna_data_dir", raising=False)
#     monkeypatch.delenv("pydna_log_dir", raising=False)

#     monkeypatch.setattr("appdirs.user_config_dir", pydna_config_dir_name)
#     monkeypatch.setattr("appdirs.user_data_dir", pydna_data_dir_name)
#     monkeypatch.setattr("appdirs.user_log_dir", pydna_log_dir_name)

#     import pydna
#     from importlib import reload

#     reload(pydna)
#     assert os.getenv("pydna_config_dir") == pydna_base_dir
#     assert os.getenv("pydna_data_dir") == pydna_base_dir
#     assert os.getenv("pydna_config_dir") == pydna_base_dir
#     pydnaenv = pydna.get_env().get_string()

#     pydna_env_vars = [v for v in os.environ if v.startswith("pydna")]

#     for envvar in pydna_env_vars:
#         assert envvar in pydnaenv
#         assert os.environ[envvar] in pydnaenv

#     subp = mock.MagicMock()
#     monkeypatch.setattr("sys.platform", "linux")
#     monkeypatch.setattr("subprocess.run", subp)

#     pydna.open_current_folder()
#     subp.assert_called_with(["xdg-open", os.getcwd()])
#     pydna.open_cache_folder()
#     subp.assert_called_with(["xdg-open", pydna_base_dir])
#     pydna.open_config_folder()
#     subp.assert_called_with(["xdg-open", pydna_base_dir])
#     pydna.open_log_folder()
#     subp.assert_called_with(["xdg-open", pydna_base_dir])
#     monkeypatch.setattr("sys.platform", "win32")
#     pydna.open_current_folder()
#     subp.assert_called_with(["start", os.getcwd()], shell=True)
#     monkeypatch.setattr("sys.platform", "darwin")
#     pydna.open_current_folder()
#     subp.assert_called_with(["open", os.getcwd()])


# def test_read_ini_file():
#     import pydna

#     pydna


# def test_without_dependency():
#     import sys
#     from unittest.mock import patch

#     with patch.dict(sys.modules, {"PIL": None}):
#         from importlib import reload

#         reload(sys.modules["pydna"])
#         import pydna

#         assert "PIL" in pydna._missing


# def test_with_dependencies():
#     pytest.importorskip("PIL")
#     pytest.importorskip("scipy")
#     pytest.importorskip("numpy")
#     import sys
#     import pydna
#     from importlib import reload

#     reload(sys.modules["pydna"])
#     import pydna

#     assert pydna._missing == []


# def test_no_xdg_open(monkeypatch):
#     subp = mock.MagicMock(side_effect=OSError(["xdg-open", os.getcwd()]))
#     monkeypatch.setattr("sys.platform", "linux")
#     monkeypatch.setattr("subprocess.run", subp)

#     import pydna

#     pydna.open_current_folder()
#     subp.assert_called_with(["xdg-open", os.getcwd()])


if __name__ == "__main__":
    pytest.main([__file__, "-vv", "-s"])
