import pandas
import pytest
from loguru import logger
import random
from muller import dataio, widgets
from muller.clustering.genotype_reorder import SortGenotypeTableWorkflow

from tests import filenames


@pytest.fixture
def sorter() -> SortGenotypeTableWorkflow:
	return SortGenotypeTableWorkflow(
		dlimit = .03,
		flimit = .97,
	)


def helper_for_testing_tables(sorter, filename):
	# Since we removed the "unsorted" genotype table, generate a random version instead.
	expected_table = dataio.import_table(filename, sheet_name = "genotype", index = "Genotype")
	unsorted_table = dataio.import_table(filename, sheet_name = "unsorted", index = "Genotype")
	expected_table = expected_table[widgets.get_numeric_columns(expected_table.columns)]
	unsorted_table = unsorted_table[widgets.get_numeric_columns(unsorted_table.columns)]
	unsorted_table.columns = [int(i) for i in unsorted_table.columns ]
	expected_table.columns = [int(i) for i in expected_table.columns]
	result = sorter.run(unsorted_table)
	result.index.name = 'Genotype'

	return result, expected_table


@pytest.mark.parametrize("filename", sorted(filenames.generic_tables.values()))
def test_generic_tables_are_sorted_correctly(sorter, filename):

	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)


@pytest.mark.parametrize("filename", sorted(filenames.model_tables.values()))
def test_model_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)


@pytest.mark.parametrize("filename", [filenames.real_tables['nature12344']])
def test_real_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)
