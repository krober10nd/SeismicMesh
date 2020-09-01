import gzip
import pathlib
import shutil

import exdown
import pytest
import requests

this_dir = pathlib.Path(__file__).resolve().parent


@pytest.mark.serial
@pytest.mark.parametrize(
    "string,lineno",
    exdown.extract(this_dir.parent / "README.md", syntax_filter="python"),
)
def test_readme(string, lineno):

    # download
    url = "http://s3.amazonaws.com/open.source.geoscience/open_data/bpvelanal2004/vel_z6.25m_x12.5m_exact.segy.gz"
    filename = url.split("/")[-1]
    with open(filename, "wb") as f:
        r = requests.get(url)
        f.write(r.content)

    # un-compress
    with gzip.open("vel_z6.25m_x12.5m_exact.segy.gz", "r") as f_in, open(
        "vel_z6.25m_x12.5m_exact.segy", "wb"
    ) as f_out:
        shutil.copyfileobj(f_in, f_out)

    try:
        # https://stackoverflow.com/a/62851176/353337
        exec(string, {"__MODULE__": "__main__"})
    except Exception:
        print(f"README.md (line {lineno}):\n```\n{string}```")
        raise
