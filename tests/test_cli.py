"""
test_cli.py

"""
import os

from click.testing import CliRunner

from aav.cli import convert

_lumi_317_path = os.path.join("tests", "data", "lumi_317_test.txt")
_lumi_370_path = os.path.join("tests", "data", "lumi_370_test.txt")
_affy_path = os.path.join("tests", "data", "affy_test.txt")
_cytoscan_path = os.path.join("tests", "data", "cytoscan_test.txt")
_open_array_path = os.path.join("tests", "data", "open_array_test.txt")


def test_no_args():
    runner = CliRunner()
    result = runner.invoke(convert, [])
    assert result.exit_code != 0


def test_readers_grch37():
    for x in [_lumi_370_path, _lumi_317_path, _affy_path, _cytoscan_path]:
        runner = CliRunner()
        result = runner.invoke(
            convert, ["-p", x, "-b", "GRCh37", "-s", "test"]
        )
        assert result.exit_code == 0
