import unittest

from muller.widgets import *


class TestWidgets(unittest.TestCase):
	def test_generate_random_color(self):
		color = generate_random_color()
		self.assertRegex(color, "#[0-9A-F]{6}")

	def test_get_numeric_columns(self):
		self.assertListEqual(['1', '66', '0', 'X9', 'x33'], get_numeric_columns(['1', '66', '0', 'X9', 'xc', 'x33', 'col4']))

	def test_generate_genotype_palette(self):
		genotypes = ['genotype-1', 'genotype-33', 'gen-5']
		expected_palette = {
			'genotype-1':  '#e6194b',
			'genotype-33': '#ffe119',
			'gen-5':       '#3cb44b',
			'genotype-0':  '#333333',
			'removed':     '#000000'
		}

		self.assertDictEqual(expected_palette, generate_genotype_palette(genotypes))

	def test_map_trajectories_to_genotype(self):
		table = pandas.DataFrame(
			{
				'genotype': ['A', 'B', 'C'],
				'members': ['A1|A2|A3', 'B1|B2', 'C1']
			}
		)
		table = table.set_index('genotype')
		expected_map = {'A1': 'A', 'A2': 'A', 'A3': 'A', 'B1':'B', 'B2': 'B', 'C1': 'C'	}
		output = map_trajectories_to_genotype(table['members'])
		self.assertDictEqual(expected_map, output)


if __name__ == "__main__":
	unittest.main()
