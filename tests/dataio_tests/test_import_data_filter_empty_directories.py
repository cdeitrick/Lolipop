from pathlib import Path

import pandas
from muller.dataio import import_tables
from loguru import logger

DATA_FOLDER = Path(__file__).parent.parent / "data"

def test_filter_empty_trajectories():
	input_column_0 = ['genotype-1', 'genotype-2', 'genotype-3', 'genotype-4', 'genotype-5', 'genotype-6']
	input_column_1 = [0.000, 0.000, 0.000, 0.111, 0.000, 0.000]
	input_column_2 = [0.000, 0.380, 0.000, 0.222, 0.000, 0.000]
	input_column_3 = [0.261, 0.432, 0.000, 0.333, 0.000, 0.000]
	input_column_4 = [1.000, 0.432, 0.000, 0.444, 1.470, 0.272]

	expected_column_0 = ['genotype-1', 'genotype-2', 'genotype-4', 'genotype-5', 'genotype-6']
	expected_column_1 = [0.000, 0.000, 0.111, 0.000, 0.000]
	expected_column_2 = [0.000, 0.380, 0.222, 0.000, 0.000]
	expected_column_3 = [0.261, 0.432, 0.333, 0.000, 0.000]
	expected_column_4 = [1.000, 0.432, 0.444, 1.470, 0.272]


	# Convert to a dataframe
	input_dataframe_definition = {
		'Genotype': input_column_0,
		0:input_column_1,
		1:input_column_2,
		2:input_column_3,
		3:input_column_4,
	}

	expected_dataframe_definition = {
		'Genotype': expected_column_0,
		0: expected_column_1,
		1: expected_column_2,
		2: expected_column_3,
		3: expected_column_4
	}
	logger.debug(input_dataframe_definition)
	input_dataframe = pandas.DataFrame(input_dataframe_definition)
	input_dataframe = input_dataframe.set_index('Genotype')
	logger.debug(input_dataframe)
	expected_dataframe = pandas.DataFrame(expected_dataframe_definition).set_index('Genotype')

	logger.debug(input_dataframe.to_string())
	result = import_tables.filter_empty_trajectories(input_dataframe)
	logger.debug(result.to_string())

	assert list(result.columns) == list(expected_dataframe.columns)

	assert len(result) == len(expected_dataframe)

	assert list(result.index) == list(expected_dataframe.index)


	#pandas.testing.assert_frame_equal(result, expected_dataframe)