import unittest

from muller.widgets import *


class TestWidgets(unittest.TestCase):

	def test_get_numeric_columns(self):
		self.assertListEqual(['1', '66', '0', 'X9', 'x33'], get_numeric_columns(['1', '66', '0', 'X9', 'xc', 'x33', 'col4']))

	def test_map_trajectories_to_genotype(self):
		table = pandas.DataFrame(
			{
				'genotype': ['A', 'B', 'C'],
				'members':  ['A1|A2|A3', 'B1|B2', 'C1']
			}
		)
		table = table.set_index('genotype')
		expected_map = {'A1': 'A', 'A2': 'A', 'A3': 'A', 'B1': 'B', 'B2': 'B', 'C1': 'C'}
		output = map_trajectories_to_genotype(table['members'])
		self.assertDictEqual(expected_map, output)


def test_get_detected_points():
	left = pandas.Series([0, 1, 1, 4, 5])
	right = pandas.Series([.23, .14, .13, 0, 0])
	result = get_detected_points(left, right, 0.03)
	assert [0, 1, 2, 3, 4] == list(result.index)

	left = pandas.Series([0, 1, 0, 0.2, 0])
	right = pandas.Series([0, .14, 0, 0, 0])
	result = get_detected_points(left, right, 0.03)
	assert [1, 2, 3] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, .23, 0, 0])
	result = get_detected_points(left, right, 0.03)
	assert [1, 2] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 0, 0])
	result = get_detected_points(left, right, 0.03)
	assert [1] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 1, 1])
	result = get_detected_points(left, right, 0.03)
	assert [1, 2, 3, 4] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 1, 1])
	result = get_detected_points(left, right, 0.03, 0.97)
	assert [1] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0, 1, 1])
	right = pandas.Series([0, 0, 0, .14, .53, 1, 1])
	result = get_detected_points(left, right, 0.03, 0.97)
	assert [3, 4] == list(result.index)

	# Check the `inner` option.
	left = pandas.Series([0, 0, .3, .4, .4, .4, 1, 1])
	right = pandas.Series([0, 0, 0, .1, .1, .1, .2, 1])
	assert [2, 3, 4, 5, 6, 7] == list(get_detected_points(left, right, .03, inner = False).index)
	assert [2, 3, 4, 5, 6] == list(get_detected_points(left, right, .03, .97, inner = False).index)
	assert [3, 4, 5, 6, 7] == list(get_detected_points(left, right, .03, inner = True).index)
	assert [3, 4, 5] == list(get_detected_points(left, right, .03, .97, inner = True).index)


if __name__ == "__main__":
	unittest.main()
