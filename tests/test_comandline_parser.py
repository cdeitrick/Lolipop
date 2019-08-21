""" A suite of tests just to see if any unexpected arguments cause the scripts to fail in weird ways."""
from pathlib import Path

import pytest

from muller import muller_workflow
from muller.commandline_parser import create_parser

DATA_FOLDER = Path(__file__).parent / "data"


@pytest.fixture
def output_folder(tmpdir) -> Path:
	return tmpdir / "output_folder"


@pytest.mark.skip
def test_minimum_commandline_parser(output_folder):
	filename = DATA_FOLDER / "tables_input_trajectories" / "5_genotypes.timeseries.tsv"
	args = ["--input", str(filename), "--output", str(output_folder)]
	args = create_parser().parse_args(args)
	w = muller_workflow.MullerWorkflow(args)
	w.run(args.filename, args.output_folder)


@pytest.mark.parametrize("cutoff", [.11, .25, .77])
def test_similiarity_cutoff(output_folder, cutoff):
	filename = DATA_FOLDER / "tables_input_trajectories" / "5_genotypes.timeseries.tsv"
	args = ["lineage", "--input", str(filename), "--output", str(output_folder), '--similarity-cutoff', cutoff]
	args = create_parser().parse_args([str(i) for i in args])
	muller_workflow.MullerWorkflow(args).run(args.filename, args.output_folder)


@pytest.mark.skip
def test_twostep_method(output_folder):
	filename = DATA_FOLDER / "tables_input_trajectories" / "5_genotypes.timeseries.tsv"
	args = ["--input", str(filename), "--output", str(output_folder), "--method", 'twostep']
	args = create_parser().parse_args([str(i) for i in args])
	mw = muller_workflow.MullerWorkflow(args)
	mw.run(args.filename, args.output_folder)
