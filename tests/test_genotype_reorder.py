from pathlib import Path

import pandas
import pytest
from loguru import logger

from muller import dataio
from muller.inheritance import SortGenotypeTableWorkflow


@pytest.fixture
def sorter() -> SortGenotypeTableWorkflow:
	return SortGenotypeTableWorkflow(.03, .15, .97, breakpoints = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0])


DATA_FOLDER = Path(__file__).parent / "data" / "tables"

generic_tables = [
	DATA_FOLDER / "generic.coexistinglineages.xlsx",
	DATA_FOLDER / "generic.genotypes.3.xlsx",
	DATA_FOLDER / "generic.genotypes.5.xlsx",
	DATA_FOLDER / "generic.genotypes.10.xlsx",
	DATA_FOLDER / "generic.small.xlsx"
]

model_tables = [
	DATA_FOLDER / "model.clonalinterferance.xlsx",
	DATA_FOLDER / "model.periodicselection.xlsx",
	DATA_FOLDER / "model.strongselection.xlsx"
]

real_tables = [
	DATA_FOLDER / "real.nature12344-s2.BYB1-G07.xlsx"
]


def helper_for_testing_tables(sorter, filename):
	original_table = dataio.import_table(filename, sheet_name = 'genotype', index = 'Genotype')
	expected_table = dataio.import_table(filename, sheet_name = 'sorted', index = 'Genotype')
	original_table.columns = [int(i) for i in original_table.columns]
	expected_table.columns = [int(i) for i in expected_table.columns]
	logger.debug(original_table.index)
	logger.debug(expected_table.index)
	result = sorter.run(original_table)
	result.index.name = 'Genotype'

	return result, expected_table


@pytest.mark.parametrize("filename", generic_tables)
def test_generic_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)


@pytest.mark.parametrize("filename", model_tables)
def test_model_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)


@pytest.mark.parametrize("filename", real_tables)
def test_real_tables_are_sorted_correctly(sorter, filename):
	result, expected_table = helper_for_testing_tables(sorter, filename)
	pandas.testing.assert_frame_equal(result, expected_table)
