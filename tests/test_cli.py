"""
test_cli.py

"""
from pathlib import Path
from click.testing import CliRunner

from aav.cli import convert

_lumi_317_path = (
    Path(__file__).parent / Path("data") / Path("lumi_317_test.txt")
)
_lumi_370_path = (
    Path(__file__).parent / Path("data") / Path("lumi_370_test.txt")
)
_affy_path = (
    Path(__file__).parent / Path("data") / Path("affy_test.txt")
)
_cytoscan_path = (
    Path(__file__).parent / Path("data") / Path("cytoscan_test.txt")
)


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
