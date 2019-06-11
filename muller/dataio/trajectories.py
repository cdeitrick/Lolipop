from pathlib import Path
from typing import Any, Optional, Tuple, Union
import math
import pandas
from loguru import logger
import re
try:
	from muller.widgets import get_numeric_columns
	from muller.dataio import import_table
except ModuleNotFoundError:
	from ..widgets import get_numeric_columns
	from . import import_table

IOTYPE = Union[str, Path]


def _convert_to_integer(value: Any, default: Optional[int] = None) -> int:
	""" Attempts to convert the input value to an integer. Returns `default` otherwise."""
	if isinstance(value, str) and (value.startswith('x') or value.startswith('X')):
		value = value[1:]
	try:
		result = int(value)
	except (TypeError, ValueError):
		result = default
	return result


def _correct_math_scale(old_data: pandas.DataFrame) -> pandas.DataFrame:
	""" Ensures the time table columns contain values between 0 and 1 and are of type `float`"""
	new_data = old_data.copy(deep = True)
	for column in old_data.columns:
		try:
			maximum_value = old_data[column].max()
		except TypeError as exception:
			logger.error(f"Some of the values in the column '{column}' could not be read as numbers.")
			logger.error(f"The column had values of {old_data[column].values}")
			raise exception
		if maximum_value > 1.0:
			logger.warning(f"The column `{column}` had values greater than 1.0. It will be converted to a float between 0 and 1.")
			new_column = old_data[column] / 100
		else:
			new_column = old_data[column].astype(float)
		new_data[column] = new_column
	return new_data
def convert_string_to_number(value:str)->float:

	pattern = "[\d]+(?:[.][\d]+)?"
	if isinstance(value, str):
		match = re.search(pattern, value)
		try:
			result = float(match.group(0))
		except AttributeError:
			# No matches were found
			logger.error(f"Could not convert '{value}' to a number.")
			result = math.nan
		except ValueError:
			# The matched string could not be converted to a number
			logger.error(f"Could not convert '{value}' to a number.")
			result = math.nan
	else:
		result = value
	return result

def _fix_column_datatypes(table: pandas.DataFrame)->pandas.DataFrame:
	# Try to extract numerical data from the table.

	for column in table.columns:
		converted_column = table[column].apply(convert_string_to_number)
		table[column] = converted_column
	return table
def _parse_table(raw_table: pandas.DataFrame, key_column: str) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	"""
		Converts column headers to integers and moves all non-integer columns to a separate dataframe.
	"""
	# Make sure the column with the series names is the index of the table.

	raw_table = raw_table.sort_values(by = key_column)
	raw_table[key_column] = [str(i) for i in raw_table[key_column].tolist()]
	raw_table.set_index(key_column, inplace = True)
	# Extract the columns which indicate timepoints of observations. Should be integers.
	frequency_columns = get_numeric_columns(raw_table.columns)

	# Extract the columns with the trajectory identifiers and frequencies at each timepoint.
	time_table = raw_table[frequency_columns]
	# Convert the columns into integers. Makes them easier to work with if they have a standard raw_table type.
	converted_frequency_columns = [_convert_to_integer(i) for i in frequency_columns]
	time_table.columns = converted_frequency_columns
	# Sort the table columns
	time_table = time_table[sorted(time_table.columns)]

	try:
		# Test to see if the table was correctly converted to numerical datatypes.
		time_table[time_table.columns[0]].max()
	except TypeError:
		logger.warning(f"Some of the values could not be interpreted as numerical data. Attempting to extract numerical data from these values...")
		time_table = _fix_column_datatypes(time_table)

	# Drop any trajectories which are never detected.
	# noinspection PyUnresolvedReferences
	time_table = time_table[(time_table.T != 0).any()]

	# Make sure the table values are between 0 and 1 and are of type `float`
	time_table = _correct_math_scale(time_table)

	# Make sure the time table contains a column for timepoint `0`.
	if 0 not in time_table.columns:
		time_table[0] = 0.0  # Should be a float to match the dtype of the other columns
		logger.warning("Warning: The input table did not have values for timepoint 0. Adding 0% for each trajectory at timepoint 0")
	# Make sure there are no missing values.
	time_table = time_table.fillna(0)
	# Extract metadata for each series.
	info_table = raw_table[[i for i in raw_table.columns if i not in frequency_columns or i == key_column]]
	return time_table, info_table


def parse_genotype_table(filename: Path, sheet_name: str = 'Sheet1') -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	""" Imports a table that lists pre-computed genotypes rather than trajectories."""
	data = import_table(filename, sheet_name = sheet_name)
	# For some reason some tables are annotated with 'genotype   ' with extra spaces.
	data.columns = [(i.strip() if isinstance(i, str) else i) for i in data.columns]
	if 'Genotype' in data.columns:
		key_column = 'Genotype'
	elif 'Unnamed: 0' in data.columns:
		key_column = 'Unnamed: 0'
	elif 'Trajectory' in data.columns:
		key_column = 'Trajectory'
	else:
		message = f"One of the columns needs to be labeled `Genotype`. Got {data.columns} instead from {filename}."
		raise ValueError(message)

	genotype_timeseries, genotype_info = _parse_table(data, key_column)
	# Make sure the genotype labels are prefixed with 'genotype-'
	if not re.search("genotype-[\d]+", genotype_timeseries.index[0]):
		genotype_timeseries.index.name = 'originalLabel'
		genotype_timeseries['Genotype'] = [f'genotype-{i}' for i in range(1, len(genotype_timeseries) + 1)]
		genotype_timeseries = genotype_timeseries.reset_index()
		genotype_timeseries = genotype_timeseries.set_index('Genotype')
		genotype_timeseries.pop('originalLabel')

	# Try to sort the genotypes by label, if posible.
	# Note: designed to sort labels of the form `genotype-[\d]+`
	try:
		sorted_index = sorted(genotype_timeseries.index, key = lambda s: float(s.split('-')[-1]))
	except ValueError:
		sorted_index = genotype_timeseries.index

	genotype_timeseries = genotype_timeseries.loc[sorted_index]

	# Remove extraneous whitespace.
	return genotype_timeseries, genotype_info


def parse_trajectory_table(filename: IOTYPE, sheet_name = 'Sheet1') -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	"""
		Reads an excel or csv file. Assumes that the file has a `Trajectory` column and a column for each timepoint.
	Parameters
	----------
	filename: Path
		The table containing the trajectories and associated metadata. Can be an excel sheet or comma/tab delimited file.
	sheet_name: str; Default 'Sheet1'
		Indicates which sheet contains the data, if an excel table is given.
	Returns
	-------
	pandas.DataFrame, pandas.DataFrame
		A timeseries dataframe
			- Index -> str
				Names unique to each trajectory.
			- Columns -> int
				The timeseries points will correspond to the frequencies for each trajectory included with the input sheet.
				Each trajectory/timepoint will include the observed frequency at each timepoint.
		A metadata dataframe
			- Index -> str
				Identical index to the timeseries dataframe.
			- Columns -> str
				All columns from the original input table that do no correspond to timepoints.
	"""

	# Read in the data table.
	raw_data = import_table(filename, sheet_name)

	key_column = 'Trajectory'
	timeseries, info = _parse_table(raw_data, key_column)

	if 'genotype' in info:
		# This file was generated by a previous run.
		info.pop('genotype')
	return timeseries, info
