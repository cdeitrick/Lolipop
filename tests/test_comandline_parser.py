""" A suite of tests just to see if any unexpected arguments cause the scripts to fail in weird ways."""
from pathlib import Path

import pytest

from muller import muller_workflow
from muller.commandline_parser import create_parser
from . import filenames
TABLE_FILENAME = filenames.generic_tables_with_trajectories['generic.genotypes.5'] # Only need to test one.

@pytest.fixture
def output_folder(tmpdir) -> Path:
	return tmpdir / "output_folder"


@pytest.mark.skip
def test_minimum_commandline_parser(output_folder):
	args = ["--input", str(TABLE_FILENAME), "--output", str(output_folder)]
	args = create_parser().parse_args(args)
	w = muller_workflow.MullerWorkflow(args)
	w.run(args.filename, args.output_folder)


@pytest.mark.parametrize("cutoff", [.11, .25, .77])
def test_similiarity_cutoff(output_folder, cutoff):
	args = ["lineage", "--input", str(TABLE_FILENAME), "--output", str(output_folder), '--pvalue', str(cutoff), "--sheetname", "trajectory"]
	args = create_parser().parse_args([str(i) for i in args])
	muller_workflow.MullerWorkflow(args).run(args.filename, args.output_folder)


@pytest.mark.skip
def test_twostep_method(output_folder):
	args = ["lineage", "--input", str(TABLE_FILENAME), "--output", str(output_folder), "--method", 'twostep', "--sheetname", "trajectory"]
	args = create_parser().parse_args([str(i) for i in args])
	mw = muller_workflow.MullerWorkflow(args)
	mw.run(args.filename, args.output_folder)
