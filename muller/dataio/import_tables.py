import io
from pathlib import Path
from typing import Optional, Union

import pandas
from loguru import logger

from muller import widgets


# noinspection PyProtectedMember
def _import_table_from_path(filename: Path, sheet_name: Optional[str] = None, index: Optional[str] = None) -> pandas.DataFrame:
	""" Imports a file as a pandas.DataFrame. Infers filetype from the filename extension/suffix.
	"""
	if filename.suffix in {'.xls', '.xlsx'}:
		data: pandas.DataFrame = pandas.read_excel(str(filename), sheet_name = sheet_name)
	else:
		sep = '\t' if filename.suffix in {'.tsv', '.tab'} else ','
		try:
			data: pandas.DataFrame = pandas.read_csv(str(filename), sep = sep)
		except UnicodeDecodeError:
			from bs4 import UnicodeDammit
			from io import StringIO
			contents = filename.read_bytes()
			unicode_contents = UnicodeDammit(contents)
			data: pandas.DataFrame = pandas.read_csv(StringIO(unicode_contents.unicode_markup), sep = sep)

	if index and index in data.columns:
		data = data.set_index(index)

	return data


def _import_table_from_string(string: str, delimiter: Optional[str] = None, index: Optional[str] = None) -> pandas.DataFrame:
	""" Imports a table represented as a basic string object."""
	# Remove unwanted whitespace.
	string = '\n'.join(i.strip() for i in string.split('\n') if i)
	if not delimiter:
		delimiter = '\t' if '\t' in string else ','
	result = pandas.read_csv(io.StringIO(string), sep = delimiter, index_col = False)
	if index:
		# Using `index_col` in `read_table()` doesn't work for some reason.
		result[index] = result[index].astype(str)
		result.set_index(index, inplace = True)
	return result

def filter_empty_trajectories(data:pandas.DataFrame)->pandas.DataFrame:
	# If the inputfiles have whitespace characters in a line they'll be imported as additional trajectories with 0% at every timepoint.
	# So basically remove any trajectories which are 0% at all timepoints.
	numeric_columns = widgets.get_numeric_columns(data.columns)
	index_to_keep = list()
	for index, row in data.iterrows():
		numeric = row[numeric_columns]
		if numeric.sum() > 0:
			index_to_keep.append(index)
	data = data.loc[index_to_keep]

	return data

def import_table(input_table: Union[str, Path], sheet_name: Optional[str] = None, index: Optional[str] = None) -> pandas.DataFrame:
	if isinstance(input_table, Path):
		data = _import_table_from_path(input_table, sheet_name, index)
	else:
		data = _import_table_from_string(input_table, index = index)
	# Make sure the x-values are numeric


	def _cast_to_int(value)->int:
		try:
			return int(value)
		except (ValueError, TypeError):
			return value

	try:
		# Using float causes problems
		data.columns = [_cast_to_int(i) for i in data.columns]
	except ValueError:
		# The columns could not be converted to str
		logger.warning(f"The columns could not be converted to int: {data.columns}")
	try:
		# Keep str columns to the left
		sorted_columns = sorted(data.columns, key = lambda s: s if isinstance(s, int) else -1)
		data = data[sorted_columns]
	except TypeError as exception:
		logger.error(data.columns)
		raise exception



	return data
