import unittest
from pathlib import Path

import pandas
import pytest

from muller.dataio import import_table, parse_genotype_table, parse_trajectory_table
from muller.dataio.import_trajectories import _convert_to_integer, _correct_math_scale
from muller.widgets import get_numeric_columns

DATA_FOLDER = Path(__file__).parent / "data"


@pytest.fixture
def truth_trajectory_table():
	trajectory_table = """Trajectory	0	17	25	44	66	75	90
	1	0	0	0.261	1	1	1	1
	2	0	0	0	0.525	0.454	0.911	0.91
	3	0	0	0	0.147	0.45	0.924	0.887
	4	0	0	0	0	0.211	0.811	0.813
	5	0	0	0	0.403	0.489	0.057	0.08
	6	0	0	0	0	0	1	1
	7	0	0	0	0.273	0.781	1	1
	8	0	0	0	0	0.345	0.833	0.793
	9	0	0	0	0	0	0.269	0.34
	10	0	0	0.117	0	0	0	0.103
	11	0	0	0	0.108	0.151	0	0
	12	0	0.125	0	0.153	0.181	0.175	0.191
	13	0	0	0	0	0.258	0.057	0.075
	14	0	0.38	0.432	0	0	0	0
	15	0	0	0.066	0.104	0.062	0	0
	16	0	0	0	0	0.209	0.209	0
	17	0	0	0	0	0	0.266	0.312
	18	0	0	0	0.115	0	0.131	0
	19	0	0	0	0.188	0.171	0.232	0.244
	20	0	0	0	0.138	0.295	0	0.081
	21	0	0	0	0.114	0	0.11	0.123
	"""
	table = import_table(trajectory_table, index = 'Trajectory')
	table.columns = [int(i) for i in table.columns]
	table = table.astype(float)
	return table


@pytest.fixture
def truth_genotype_table():
	genotype_table = """Genotype	0	17	25	44	66	75	90
	genotype-1	0	0	0.261	1	1	1	1
	genotype-2	0	0.38	0.432	0	0	0	0
	genotype-3	0	0	0	0	0	1	1
	genotype-4	0	0	0	0.525	0.454	0.911	0.91
	genotype-5	0	0	0	0.147	0.45	0.924	0.887
	genotype-6	0	0	0	0.273	0.781	1	1
	genotype-7	0	0	0	0.188	0.171	0.232	0.244
	genotype-8	0	0	0	0.403	0.489	0.057	0.08
	genotype-9	0	0	0.117	0	0	0	0.103
	genotype-10	0	0	0	0.138	0.295	0	0.081
	genotype-11	0	0	0	0	0.278	0.822	0.803
	genotype-12	0	0	0	0	0.2335	0.133	0.0375
	genotype-13	0	0	0.033	0.106	0.1065	0	0
	genotype-14	0	0	0	0	0	0.2675	0.326
	genotype-15	0	0	0	0.1145	0	0.1205	0.0615
	"""
	table = import_table(genotype_table, index = 'Genotype')
	table.columns = [int(i) for i in table.columns]
	table = table.astype(float)
	return table


@pytest.mark.parametrize('filename', list((DATA_FOLDER / "test_trajectory_tables").iterdir()))
def test_import_trajectory_table(filename, truth_trajectory_table):
	test_table, info = parse_trajectory_table(filename)
	pandas.testing.assert_frame_equal(truth_trajectory_table, test_table)


@pytest.mark.parametrize(
	'filename',
	[
		"B1_muller_try1.muller_genotypes.original.tsv",
		"B1_muller_try1.muller_genotypes.original.xls",
		"B1_muller_try1.muller_genotypes.original.xlsx"
	]
)
def test_import_genotype_tables(filename, truth_genotype_table):
	folder = DATA_FOLDER / "test_genotype_tables"
	test_table, info = parse_genotype_table(folder / filename)
	pandas.testing.assert_frame_equal(truth_genotype_table, test_table)


class TestImportTable(unittest.TestCase):
	def test_get_numeric_columns(self):
		columns = ["abc", '123', 456, 'X66', 'x0']
		output = get_numeric_columns(columns)

		self.assertListEqual(['123', 456, 'X66', 'x0'], output)

	def test_correct_math_scale(self):
		string = """Trajectory	X0	X1	X2	X3	X4	X5
			trajectory-A2	0	0	0	6	35	4
			trajectory-A3	0	0	0	0	45	5"""
		df = import_table(string, index = 'Trajectory')
		self.assertListEqual([35, 45], df['X4'].tolist())

		fdf = _correct_math_scale(df)
		self.assertListEqual([.35, .45], fdf['X4'].tolist())
		self.assertEqual(str(fdf['X0'].dtype), 'float64')

	def test_convert_to_integer(self):
		self.assertEqual(66, _convert_to_integer('66'))
		self.assertEqual(55, _convert_to_integer('x55'))
		self.assertEqual(123, _convert_to_integer(123))
		self.assertEqual(45, _convert_to_integer('X45'))
		self.assertIsNone(_convert_to_integer('abc123'))
		self.assertEqual(1, _convert_to_integer('abc', default = 1))

	def test_import_genotype_table_does_not_crash(self):
		df, info = parse_genotype_table(DATA_FOLDER / "3_genotypes.genotypes.tsv")


if __name__ == "__main__":
	pass
