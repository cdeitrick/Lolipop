from typing import Dict, Tuple

import pandas
import pytest

from muller.inheritance import genotype_lineage
from tests import filenames


@pytest.fixture
def lineage_workflow() -> genotype_lineage.LineageWorkflow:
	return genotype_lineage.LineageWorkflow(0.03, 0.97, 0.05)


def helper_read_table(filename: str) -> Tuple[pandas.DataFrame, Dict[str, str]]:
	""" Reads in a table that contains all relevant information pertaining to a set of genotypes."""
	table_filename = filename

	genotypes = pandas.read_excel(table_filename, sheet_name = "genotype").set_index('Genotype')
	edges = pandas.read_excel(table_filename, sheet_name = "edges")
	# Convert the `edges` dataframe into an ordinary dictionary
	edges = edges.set_index("Identity")["Parent"].to_dict()
	return genotypes, edges


@pytest.mark.parametrize(
	"filename", filenames.generic_tables_with_trajectories.values()
)
def test_lineage_of_generic_tables(lineage_workflow, filename: str):
	""" A parametrized function to assess how well the scripts infer genotypes from generic data."""

	# Note: the table with 10 genotypes is a genotype table rather than a trajectory table.
	# Only difference is what the tables will be used for, but the distinction matters when every
	# series is a genotype rather than a trajectory, since genotypes can't be clustered into a genotype.
	genotypes, expected_lineage = helper_read_table(filename)

	if len(genotypes) == 10:
		return None

	result_lineage = lineage_workflow.run(genotypes)
	result = result_lineage.clusters.as_dict()

	# A helper function for debugging dictionary comparisons.
	filenames.compare_dictionaries(result, expected_lineage)
	assert result == expected_lineage
