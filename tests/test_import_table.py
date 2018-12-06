import unittest
from pathlib import Path
from muller.import_table import *
from io import StringIO
FILENAME = Path(__file__).parent / "data" / "B1_muller_try1.xlsx"
test_table = pandas.read_excel(FILENAME)
table_text = """0,17,25,44,66,75,90
0,0.0,0.261,1.0,1.0,1.0,1.0
0,0.0,0.0,0.525,0.45399999999999996,0.9109999999999999,0.91
0,0.0,0.0,0.147,0.45,0.924,0.887
0,0.0,0.0,0.0,0.21100000000000002,0.8109999999999999,0.813
0,0.0,0.0,0.40299999999999997,0.489,0.057,0.08
0,0.0,0.0,0.0,0.0,1.0,1.0
0,0.0,0.0,0.273,0.7809999999999999,1.0,1.0
0,0.0,0.0,0.0,0.345,0.833,0.7929999999999999
0,0.0,0.0,0.0,0.0,0.26899999999999996,0.34
0,0.0,0.11699999999999999,0.0,0.0,0.0,0.10300000000000001
0,0.0,0.0,0.10800000000000001,0.151,0.0,0.0
0,0.125,0.0,0.153,0.18100000000000002,0.175,0.191
0,0.0,0.0,0.0,0.258,0.057,0.075
0,0.38,0.43200000000000005,0.0,0.0,0.0,0.0
0,0.0,0.066,0.10400000000000001,0.062,0.0,0.0
0,0.0,0.0,0.0,0.209,0.209,0.0
0,0.0,0.0,0.0,0.0,0.266,0.312
0,0.0,0.0,0.115,0.0,0.131,0.0
0,0.0,0.0,0.188,0.171,0.23199999999999998,0.244
0,0.0,0.0,0.138,0.295,0.0,0.081
0,0.0,0.0,0.114,0.0,0.11,0.12300000000000001"""
truth_table = pandas.read_csv(StringIO(table_text))
truth_table.columns = [int(i) for i in truth_table.columns]
class TestImportTable(unittest.TestCase):
	def test_get_numeric_columns(self):
		columns = ["abc", '123', 456, 'X66']
		output = get_numeric_columns(columns)

		self.assertListEqual(['123', 456, 'X66'], output)

	def test_correct_math_scale(self):
		output = correct_math_scale(test_table)
		output = output[[0,17,25,44,66,75,90]]
		output.columns = truth_table.columns
		pandas.testing.assert_frame_equal(truth_table, output)

	def test_import_table(self):
		output, info = import_trajectory_table(FILENAME)
		pandas.testing.assert_frame_equal(truth_table, output)
if __name__ == "__main__":
	pass