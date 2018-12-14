import unittest
from import_table import import_table_from_string
import math
from muller_genotypes.similarity import *
class TestSimilarity(unittest.TestCase):
	def test_calculate_p_value(self):
		string = """
		Trajectory	0	1	2	3	4	5
		trajectory-A1	0	0	0	0.1	0.5	0.5
		trajectory-B1	0	0.1	0.15	0.03	0	0
		trajectory-A2	0	0	0	0.06	0.35	0.4
		trajectory-C1	0	0	0	0.3	0.7	1
		trajectory-A3	0	0	0	0	0.45	0.5
		trajectory-B2	0	0.07	0.1	0.02	0.01	0
		trajectory-F1	0	0	1	1	1	1
		trajectory-F2	0	1	0	0	0	0
		trajectory-F3	0	1	1	1	1	1
		trajectory-U1	0	0	0	0	0	0
		trajectory-U2	0	0.02	0.01	0.01	0.02	0
		"""
		table = import_table_from_string(string, index = 'Trajectory')

		# Check when both are undetected.
		p_value = calculate_p_value(table.loc['trajectory-U1'], table.loc['trajectory-U2'], .05, .95)
		self.assertEqual(1.0, p_value.pvalue)
		self.assertIsInstance(p_value.pvalue, float)

		# check when both are fixed or undetected, but do not overlap.
		p_value = calculate_p_value(table.loc['trajectory-F1'], table.loc['trajectory-F2'], .05, .95)
		self.assertEqual(0.0, p_value.pvalue)

		# Check when both are undetected but have overlap.
		calculation = calculate_p_value(table.loc['trajectory-F1'], table.loc['trajectory-F3'], .05, .95)
		self.assertEqual(1.0, calculation.pvalue)

	def test_calculate_pairwise_similarity(self):
		pass
if __name__ == "__main__":
	unittest.main()