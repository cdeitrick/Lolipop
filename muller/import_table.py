from pathlib import Path
from typing import Any, Optional, Tuple

import pandas

from widgets import get_numeric_columns


def correct_math_scale(old_data: pandas.DataFrame) -> pandas.DataFrame:
	""" Ensures the time table columns contain values between 0 and 1 and are of type `float`"""
	new_data = old_data.copy(deep = True)
	for column in old_data.columns:
		if max(old_data[column]) > 1.0:
			new_column = old_data[column] / 100
		else:
			new_column = old_data[column].astype(float)
		new_data[column] = new_column
	return new_data


def import_table(filename: Path, sheet_name: str) -> pandas.DataFrame:
	""" Imports a file as a pandas.DataFrame. Infers filetype from the filename extension/suffix.
	"""
	if filename.suffix in {'.xls', '.xlsx'}:
		data: pandas.DataFrame = pandas.read_excel(str(filename), sheet_name = sheet_name)

	else:
		sep = '\t' if filename.suffix in {'.tsv', '.tab'} else ','
		data: pandas.DataFrame = pandas.read_table(str(filename), sep = sep)

	data = data[sorted(data.columns, key = lambda s: str(s))]
	return data


def import_genotype_table(filename: Path, sheetname: str) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	""" Imports a table that lists pre-computed genotypes rather than trajectories."""
	data = import_table(filename, sheet_name = sheetname)

	if 'Genotype' in data.columns:
		key_column = 'Genotype'
	elif 'Unnamed: 0' in data.columns:
		key_column = 'Unnamed: 0'
	else:
		message = "One of the columns needs to be labeled `Genotype`"
		raise ValueError(message)

	if 'members' not in data.columns:
		message = "The genotype must have a 'members' column with the names of all trajectories contained in the genotype. Individual trajectory names must be separated by '|'"
		raise ValueError(message)

	genotype_timeseries, genotype_info = _parse_table(data, key_column)
	return genotype_timeseries, genotype_info


def _convert_to_integer(value: Any, default: Optional[int] = None) -> int:
	""" Attempts to convert the input value to an integer. Returns `default` otherwise."""
	if isinstance(value, str) and (value.startswith('x') or value.startswith('X')):
		value = value[1:]
	try:
		result = int(value)
	except (TypeError, ValueError):
		result = default
	return result


def _parse_table(raw_table: pandas.DataFrame, key_column: str) -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	"""
		Converts column headers to integers and moves all non-integer columns to a separate dataframe.
	"""
	# Make sure the column with the series names is the index of the table.
	raw_table[key_column] = [str(i) for i in raw_table[key_column].tolist()]
	raw_table.set_index(key_column, inplace = True)

	# Extract the columns which indicate timepoints of observations. Should be integers.
	frequency_columns = get_numeric_columns(raw_table.columns)

	# Extract the columns with the trajectory identifiers and frequencies at each timepoint.
	time_table = raw_table[frequency_columns]
	# Convert the columns into integers. Makes them easier to work with if they have a standard raw_table type.
	converted_frequency_columns = [_convert_to_integer(i) for i in frequency_columns]
	time_table.columns = converted_frequency_columns

	# Drop any trajectories which are never detected.
	time_table = time_table[(time_table.T != 0).any()]

	# Make sure the table values are between 0 and 1 and are of type `float`
	time_table = correct_math_scale(time_table)

	# Make sure the time table contains a column for timepoint `0`.
	if 0 not in time_table.columns:
		time_table[0] = 0.0  # Should be a float to match the dtype of the other columns
		print("Warning: The input table did not have values for timepoint 0. Adding 0% for each trajectory at timepoint 0")

	# Extract metadata for each series.
	info_table = raw_table[[i for i in raw_table.columns if i not in frequency_columns]]
	return time_table, info_table


def import_trajectory_table(filename: Path, sheet_name = 'Sheet1') -> Tuple[pandas.DataFrame, pandas.DataFrame]:
	"""
		Reads an excel or csv file. Assumes that the file has the following columns:
		- Population:
		- Trajectory:
		- Position:
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
	return timeseries, info


if __name__ == "__main__":
	path = Path(__file__).parent.parent / "tests" / "data" / "2_genotypes_with_string_columns.tsv"
	data, _ = import_trajectory_table(path)
	print(data.to_string())
