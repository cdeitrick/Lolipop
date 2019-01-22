import unittest

from muller_genotypes.filters import *
from muller_genotypes.filters import _get_backgrounds_present_at_multiple_timepoints
from import_data import import_table_from_string


class TestGenotypeFilters(unittest.TestCase):
	def test_get_fuzzy_backgrounds(self):
		string_table = """	
			Genotype	0	1	2	3	4	5
			Trajectory-A	0	0	0	0.1	0.5	0.5
			Trajectory-B	0	0.1	0.15	0.03	0	0
			Trajectory-C	0	0	0	0.3	0.7	.92"""
		table = import_table_from_string(string_table, index = 'Genotype')

		cutoffs = [1, .9, .8]
		expected_background = [0, 0, 0, 0.3, 0.7, 0.92]
		output_background, (fdc, ffc) = get_fuzzy_backgrounds(table, cutoffs)
		self.assertEqual(1 - .9, fdc)
		self.assertEqual(.9, ffc)
		self.assertEqual(1, len(output_background))
		self.assertListEqual(expected_background, output_background.iloc[0].tolist())

		cutoffs = [1.0, .99, .95]
		self.assertRaises(ValueError, get_fuzzy_backgrounds, table, cutoffs)

	def test_get_first_timepoint(self):
		series = pandas.Series([0, 0, 0.1, 0.3, 0.7, 0.92], index = [0, 11, 65, 400, 401, 402])
		self.assertEqual(400, get_first_timpoint_above_cutoff(series, .15))
		self.assertEqual(401, get_first_timpoint_above_cutoff(series, .5))

	def test_check_if_genotype_is_invalid(self):
		background_detected_point = 7
		background_fixed_point = 12
		index = [0, 3, 7, 10, 12, 15, 20, 25]
		valid_genotype = pandas.Series([0, 0, 0, 0, 0, .1, .3, .7], index = index)

		self.assertFalse(check_if_genotype_is_invalid(valid_genotype, background_detected_point, background_fixed_point, 0.05))

		invalid_genotype = pandas.Series([0, .1, 0.2, .3, .1, .1, .1, .1], index = index)
		self.assertTrue(check_if_genotype_is_invalid(invalid_genotype, background_detected_point, background_fixed_point, 0.05))

		valid_genotype = pandas.Series([0, .1, .0, .1, .1, .1, .1, .1])
		self.assertFalse(check_if_genotype_is_invalid(valid_genotype, background_detected_point, background_fixed_point, .05))

	def test_get_invalid_genotype(self):
		string = """
		Genotype	0	1	2	3	4	5	members
		Trajectory-A	0	0	0	0.1	0.5	0.5	A1|A2|A3
		Trajectory-B	0	0.1	0.15	0.03	0	0	B1|B2
		Trajectory-C	0	0	0	0.3	0.97	1	C1
		Trajectory-D	0	.1	.1	.1	.1	.1	D1
		"""
		table = import_table_from_string(string, index = 'Genotype')
		table.pop('members')
		output = find_first_invalid_genotype(table, .05, [1, .9, .8, .7, .6])
		self.assertEqual('Trajectory-D', output)

	def test_get_backgrounds_present_at_multiple_timepoints(self):
		backgrounds = pandas.DataFrame(
			{
				'genotype': ['genotype-A', 'genotype-B', 'genotype-C'],
				0:          [0, 0.01, 0],
				2:          [0, 1, .1],
				5:          [0, 0, .5],
				7:          [.4, 0.01, .97],
				8:          [1, 0, 1],
				9:          [1, 0, 1]
			}
		).set_index('genotype')
		output = _get_backgrounds_present_at_multiple_timepoints(backgrounds, 0.05)

		self.assertListEqual(['genotype-A', 'genotype-C'], list(output.index))

if __name__ == "__main__":
	unittest.main()
