import unittest
from pathlib import Path
from muller.import_table import *
class TestImportTable(unittest.TestCase):
	def test_get_numeric_columns(self):
		columns = ["abc", '123', 456, 'X66']
		output = get_numeric_columns(columns)

		self.assertListEqual(['123', 456, 'X66'], output)

if __name__ == "__main__":
	pass