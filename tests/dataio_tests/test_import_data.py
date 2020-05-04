from pathlib import Path

import pandas
import pytest
from typing import *
from muller.dataio import import_table, parse_genotype_table, parse_trajectory_table
from muller.dataio.import_timeseries import _convert_to_integer, _correct_math_scale
from muller import widgets
from tests import filenames, twidgets
from loguru import logger

DATA_FOLDER = Path(__file__).parent.parent / "data"

def make_table(parent_folder:Path, name:str):
	filename = filenames.real_tables['B1']
	table = pandas.read_excel(filename, sheet_name = name)

	table.to_csv(parent_folder / "B1.csv", index = False)
	table.to_csv(parent_folder / "B1.tsv", sep = "\t", index = False)
	table.to_excel(parent_folder / "B1.xlsx", index = False)

	table.loc[:, name.capitalize()] = table[name.capitalize()].astype(str)
	table = table.set_index(name.capitalize())
	result = {
		'extensions': ['.csv', '.tsv', '.xlsx'],
		'basename':   parent_folder / "B1",
		'truthset':   table
	}
	return result

@pytest.fixture
def truth_trajectory_tables(tmp_path) -> Dict[str, Any]:
	result = make_table(tmp_path, 'trajectory')
	return result


@pytest.fixture
def truth_genotype_tables(tmp_path):
	result = make_table(tmp_path, 'genotype')
	return result

def _run(truth_tables, name):

	extensions = truth_tables['extensions']
	basename = truth_tables['basename']
	truthset = truth_tables['truthset']
	for extension in extensions:
		filename = basename.with_suffix(extension)
		test_table, info = parse_trajectory_table(filename)
		test_table.index.name = name.capitalize
		logger.debug(test_table.index)
		twidgets.assert_frame_equal(test_table, truthset)

def test_import_trajectory_table(truth_trajectory_tables):
	_run(truth_trajectory_tables, 'trajectory')


def test_import_genotype_tables(truth_genotype_tables):
	_run(truth_genotype_tables, 'genotype')


def test_get_numeric_columns():
	columns = ["abc", '123', 456, 'X66', 'x0']
	output = widgets.get_numeric_columns(columns)

	assert output == ['123', 456, 'X66', 'x0']


def test_correct_math_scale():
	string = """Trajectory	X0	X1	X2	X3	X4	X5
		trajectory-A2	0	0	0	6	35	4
		trajectory-A3	0	0	0	0	45	5"""
	df = import_table(string, index = 'Trajectory')
	assert df['X4'].tolist() == [35, 45]

	fdf = _correct_math_scale(df)
	assert fdf['X4'].tolist() == [0.35, 0.45]
	assert fdf['X4'].dtype == float


@pytest.mark.parametrize("value, expected", [
	('66', 66),
	('x55', 55),
	(123, 123),
	('X45', 45),
	('abc123', None)
])
def test_convert_to_integer(value, expected):
	assert _convert_to_integer(value) == expected


if __name__ == "__main__":
	pass
