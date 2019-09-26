import pandas
import pytest
from loguru import logger

from muller import dataio
from muller.inheritance import SortGenotypeTableWorkflow
from .filenames import generic_tables_with_trajectories, model_tables, real_tables


@pytest.fixture
def sorter() -> SortGenotypeTableWorkflow:
	return SortGenotypeTableWorkflow(.03, .15, .97, breakpoints = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0])


def helper_for_testing_tables(sorter, filename):
	original_table = dataio.import_table(filename, sheet_name = 'unsorted', index = 'Genotype')
	expected_table = dataio.import_table(filename, sheet_name = 'genotype', index = 'Genotype')
	original_table.columns = [int(i) for i in original_table.columns]
	expected_table.columns = [int(i) for i in expected_table.columns]
	logger.debug(original_table.index)
	logger.debug(expected_table.index)
	result = sorter.run(original_table)
	result.index.name = 'Genotype'

	return result, expected_table


@pytest.mark.parametrize("filename", sorted(generic_tables_with_trajectories.values()))
def test_generic_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)


@pytest.mark.parametrize("filename", sorted(model_tables.values()))
def test_model_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)


@pytest.mark.parametrize("filename", sorted(real_tables.values()))
def test_real_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)
