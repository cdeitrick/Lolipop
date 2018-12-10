import unittest

from muller.import_table import *
from muller.import_table import _convert_to_integer, _correct_math_scale

DATA_FOLDER = Path(__file__).parent / "data"
trajectory_table = """Trajectory	X0	X1	X2	X3	X4	X5
trajectory-A1	0	0	0	0.1	0.5	0.5
trajectory-B1	0	0.1	0.15	0.03	0	0
trajectory-A2	0	0	0	0.06	0.35	0.4
trajectory-A3	0	0	0	0	0.45	0.5
trajectory-B2	0	0.07	0.1	0.02	0.01	0"""


class TestImportTable(unittest.TestCase):
	def test_get_numeric_columns(self):
		columns = ["abc", '123', 456, 'X66', 'x0']
		output = get_numeric_columns(columns)

		self.assertListEqual(['123', 456, 'X66', 'x0'], output)

	def test_import_table_from_string(self):
		result = import_table_from_string(trajectory_table, index = 'Trajectory')
		# Keep the test simple for now. just test that the required columns and values are present.
		truth_index = [
			'trajectory-A1', 'trajectory-B1', 'trajectory-A2', 'trajectory-A3', 'trajectory-B2'
		]
		self.assertListEqual(truth_index, list(result.index))
		self.assertListEqual([0.5, 0.0, 0.4, 0.5, 0.0], result['X5'].tolist())

	def test_import_trajectory_table_from_string(self):
		result, result_info = import_trajectory_table(trajectory_table)
		self.assertTrue(result_info.empty)
		truth_index = [
			'trajectory-A1', 'trajectory-B1', 'trajectory-A2', 'trajectory-A3', 'trajectory-B2'
		]
		self.assertListEqual(truth_index, list(result.index))
		self.assertListEqual([0.5, 0.0, 0.4, 0.5, 0.0], result[5].tolist())

	def test_correct_math_scale(self):
		string = """Trajectory	X0	X1	X2	X3	X4	X5
			trajectory-A2	0	0	0	6	35	4
			trajectory-A3	0	0	0	0	45	5"""
		df = import_table_from_string(string, index = 'Trajectory')
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
		df, info = import_genotype_table(DATA_FOLDER / "3_genotypes.genotypes.tsv")


if __name__ == "__main__":
	pass
