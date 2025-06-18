#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


pytest.importorskip("requests")


@pytest.fixture
def requests_mock(request):
    import requests_mock as rm_module

    m = rm_module.Mocker()
    m.start()
    request.addfinalizer(m.stop)
    return m


def test_web(requests_mock, monkeypatch):
    from pydna import download

    # monkeypatch.setenv("pydna_cached_funcs", "")

    import io

    flo = io.BytesIO(b"some text data")

    requests_mock.get(
        "http://www.fake.com/hej.txt",
        headers={
            "last-modified": "Mon, 01 Jan 2001 00:00:00 GMT",  # 978307200
            "content-length": "14",
        },
        body=flo,
    )

    tx = download.download_text("http://www.fake.com/hej.txt")

    assert tx == "some text data"


if __name__ == "__main__":
    pytest.main([__file__, "-x", "-vv", "-s"])
