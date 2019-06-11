import io
from pathlib import Path
from typing import Optional, Union

import pandas


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


def import_table(input_table: Union[str, Path], sheet_name: Optional[str] = None, index: Optional[str] = None) -> pandas.DataFrame:
	if isinstance(input_table, Path):
		data = _import_table_from_path(input_table, sheet_name, index)
	else:
		data = _import_table_from_string(input_table, index = index)
	return data
