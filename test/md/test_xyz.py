"""Test functions handling xyz files."""
import pathlib

from dztools.md.center_xyz import center_periodic

# import numpy as np
# import pytest
# import tomli

HERE = pathlib.Path(__file__).resolve().parent

def test_center_periodic(tmp_path):
    """Check that we can modify the velocities with an engine,
    and that they are not equal to zero."""
    # folder we wil run from
    folder = tmp_path / "temp"
    folder.mkdir()

    tocenter = str(HERE) + '/files/conf.xyz'
    outfile = str(folder) + '/conf_c.xyz'

    input = ['-i', tocenter, '-o', outfile, '-c', '12.4138', '-idx', '2']
    center_periodic(input)

    rip = []
    with open(outfile) as read:
        for idx, line in enumerate(read):
            rip = line.rstrip().split()
            if idx == 3:
                break
    assert rip[0] == "C"
    for i in rip[1:]:
        assert float(i) == float(0)
