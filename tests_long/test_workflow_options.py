"""
	Runs the full workflow with various options given. There are separate tests from some of the options already (such as reading in annotations,"
	but these tests basically make sure all the options were implemented properly.
"""
from pathlib import Path
from typing import *

import pandas
import pytest

from muller import dataio, workflows


@pytest.fixture
def output(tmp_path) -> Path:
	""" Returns a path to a temporary directory which can be used as the output of each workflow."""
	return tmp_path


@pytest.fixture
def table() -> pandas.DataFrame:
	""" Returns the path to the input dataset for these tests. Marked as a fixture so that it can be called as a method parameter rather then
		as a standalone method call within the test methods.
	"""

	folder_data = Path(__file__).parent / "data"
	filename = folder_data / "generic.genotypes.10.xlsx"
	table = dataio.import_table(filename, index = 'Trajectory', sheet_name = "trajectory")

	return table


def make_known_genotypes_file(path: Path, genotypes):
	if isinstance(genotypes[0], list):
		genotypes = [",".join(i) for i in genotypes]
	filename = path / "known_genotypes.txt"
	contents = "\n".join(genotypes)
	filename.write_text(contents)

	return filename


def run_workflow(table: pandas.DataFrame, known_genotypes: Optional[Path]):
	# Implemented to make the test methods a little cleaner.
	result = workflows.run_genotype_inference_workflow(
		trajectoryio = table,
		metric = "binomial",
		dlimit = 0.03,
		slimit = 0.15,
		flimit = 0.97,
		pvalue = 0.05,
		known_genotypes = known_genotypes
	)
	return result


def test_known_genotypes_one(table, output):
	# First, run the workflow to make sure the genotypes don't already exist.

	genotype = ["genotype-green-1", "genotype-purple-1"]
	filename_known_genotypes = make_known_genotypes_file(output, [genotype])

	result = run_workflow(table, None)
	# Make sure trajectories green and purple are in separate clusters
	assert not any(set(genotype) <= set(i) for i in result.genotype_members.values())

	# Now rerun the workflow with the known genotype.
	workflowoutput = run_workflow(table, filename_known_genotypes)

	assert any(set(genotype) <= set(i) for i in workflowoutput.genotype_members.values())


def test_known_genotypes_three(table, output):
	genotypes = [
		['genotype-green-1', 'genotype-purple-1'],
		['genotype-steelblue-1', 'genotype-teal-1'],
		['genotype-orange-1', 'genotype-gold-1', 'genotype-sienna-1']
	]

	# Already made sure these genotypes didn't exist in the other method.
	# Skip running the control
	filename_known_genotypes = make_known_genotypes_file(output, genotypes)
	result = run_workflow(table, filename_known_genotypes)

	# Test whether the expected genotypes are a subset of any actual genotypes.
	result1 = [set(genotypes[0]) <= set(i) for i in result.genotype_members.values()]
	result2 = [set(genotypes[1]) <= set(i) for i in result.genotype_members.values()]
	result3 = [set(genotypes[2]) <= set(i) for i in result.genotype_members.values()]

	assert result1 and result2 and result3


def test_gene_aliases():
	pass


def test_genotype_colors():
	pass


if __name__ == "__main__":
	pass
