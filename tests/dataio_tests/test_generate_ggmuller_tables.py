import unittest

import pandas

from muller.dataio import import_table
from muller.dataio.generate_tables import add_ancestral_genotype, _compile_parent_linkage, _convert_genotype_table_to_population_table

edges_table = """
	Parent	Identity
	genotype-0	genotype-1
	genotype-1	genotype-4
	genotype-4	genotype-5
	genotype-0	genotype-2
	genotype-1	genotype-3
"""


class TestEdgesTable(unittest.TestCase):
	def test_compile_parent_linkage(self):
		edges = import_table(edges_table)
		expected_linkage = {
			'genotype-1': ['genotype-4', 'genotype-3'],
			'genotype-4': ['genotype-5'],
			'genotype-0': ['genotype-1', 'genotype-2']
		}
		output = _compile_parent_linkage(edges)
		self.assertDictEqual(expected_linkage, output)

	def test_generate_edges_table(self):
		pass


class TestPopulationTable(unittest.TestCase):
	def test_subtract_children_from_parent(self):
		pass

	def test_convert_genotype_table_to_population_table(self):
		test_table = import_table(
			"""
			Genotype	0	1	22	33
			genotype-1	0	.12	.5	0
			genotype-2	0	.23	.5	.7
			""", index = 'Genotype'
		)
		expected_table = import_table(
			"""
			Generation	Identity	Population
			0	genotype-1	0.0
			1	genotype-1	12.0
			22	genotype-1	50.0
			33	genotype-1	0.0
			0	genotype-2	0.0
			1	genotype-2	23.0
			22	genotype-2	50.0
			33	genotype-2	70.0
		"""
		)
		# expected_table['Generation'] = expected_table['Generation'].astype(str)
		test_output = _convert_genotype_table_to_population_table(test_table)
		pandas.testing.assert_frame_equal(expected_table, test_output)

	def test_append_genotype_0(self):
		basic_population_table = """
			Generation	Identity	Population
			1	genotype-a	0
			2	genotype-a	0
			3	genotype-a	0
			5	genotype-a	10
			1	genotype-b	0
			2	genotype-b	50
			3	genotype-b	55
			5	genotype-b	90
		"""
		basic_population_table = import_table(basic_population_table)
		test_population_table = add_ancestral_genotype(basic_population_table)

		expected_table = """
			Generation	Identity	Population
			1	genotype-a	0
			2	genotype-a	0
			3	genotype-a	0
			5	genotype-a	10
			1	genotype-b	0
			2	genotype-b	50
			3	genotype-b	55
			5	genotype-b	90
			1	genotype-0	100
			2	genotype-0	50
			3	genotype-0	45
			5	genotype-0	0
		"""
		expected_table = import_table(expected_table)

		pandas.testing.assert_frame_equal(expected_table, test_population_table)


if __name__ == "__main__":
	unittest.main()
