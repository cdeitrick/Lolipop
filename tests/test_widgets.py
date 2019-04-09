from unittest.mock import patch, MagicMock

import pandas
import dataio
import widgets


def test_get_numeric_columns():
	assert ['1', '66', '0', 'X9', 'x33'] == widgets.get_numeric_columns(['1', '66', '0', 'X9', 'xc', 'x33', 'col4'])


def test_map_trajectories_to_genotype():
	table = pandas.DataFrame(
		{
			'genotype': ['A', 'B', 'C'],
			'members':  ['A1|A2|A3', 'B1|B2', 'C1']
		}
	)
	table = table.set_index('genotype')
	expected_map = {'A1': 'A', 'A2': 'A', 'A3': 'A', 'B1': 'B', 'B2': 'B', 'C1': 'C'}
	output = widgets.map_trajectories_to_genotype(table['members'])
	assert expected_map == output


def test_get_detected_points():
	left = pandas.Series([0, 1, 1, 4, 5])
	right = pandas.Series([.23, .14, .13, 0, 0])
	result = widgets.get_detected_points(left, right, 0.03)
	assert [0, 1, 2, 3, 4] == list(result.index)

	left = pandas.Series([0, 1, 0, 0.2, 0])
	right = pandas.Series([0, .14, 0, 0, 0])
	result = widgets.get_detected_points(left, right, 0.03)
	assert [1, 2, 3] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, .23, 0, 0])
	result = widgets.get_detected_points(left, right, 0.03)
	assert [1, 2] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 0, 0])
	result = widgets.get_detected_points(left, right, 0.03)
	assert [1] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 1, 1])
	result = widgets.get_detected_points(left, right, 0.03)
	assert [1, 2, 3, 4] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 1, 1])
	result = widgets.get_detected_points(left, right, 0.03, 0.97)
	assert [1] == list(result.index)

	left = pandas.Series([0, 0, 0, 0, 0, 1, 1])
	right = pandas.Series([0, 0, 0, .14, .53, 1, 1])
	result = widgets.get_detected_points(left, right, 0.03, 0.97)
	assert [3, 4] == list(result.index)

	# Check the `inner` option.
	left = pandas.Series([0, 0, .3, .4, .4, .4, 1, 1])
	right = pandas.Series([0, 0, 0, .1, .1, .1, .2, 1])
	assert [2, 3, 4, 5, 6, 7] == list(widgets.get_detected_points(left, right, .03, inner = False).index)
	assert [2, 3, 4, 5, 6] == list(widgets.get_detected_points(left, right, .03, .97, inner = False).index)
	assert [3, 4, 5, 6, 7] == list(widgets.get_detected_points(left, right, .03, inner = True).index)
	assert [3, 4, 5] == list(widgets.get_detected_points(left, right, .03, .97, inner = True).index)


def test_get_valid_points_simple():
	left = pandas.Series([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 0])
	right = pandas.Series([0, 0, 0, .1, .2, .3, .3, .3, .3, 0, 0, 0])

	expected = pandas.DataFrame({
		'left':  [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1],
		'right': [0, 0, .1, .2, .3, .3, .3, .3, 0, 0],
	}, index = range(1, 11))
	result = widgets.get_valid_points(left, right, 0.03)
	pandas.testing.assert_frame_equal(expected, result)

	expected = pandas.DataFrame({
		'left':  [.1, .2, .3, .4, .5, .6, .7, .8, .9],
		'right': [0, 0, .1, .2, .3, .3, .3, .3, 0],
	}, index = range(1, 10))

	result = widgets.get_valid_points(left, right, 0.03, 0.97)
	pandas.testing.assert_frame_equal(expected, result)

	expected = pandas.DataFrame({
		'left':  [.3, .4, .5, .6, .7, .8],
		'right': [.1, .2, .3, .3, .3, .3],
	}, index = range(3, 9))
	result = widgets.get_valid_points(left, right, 0.03, 0.97, inner = True)
	pandas.testing.assert_frame_equal(expected, result)


def test_get_valid_points_complex():
	left = pandas.Series([0, 0.653, 1, 1, 1, 0.91, 0.907, 1])
	right = pandas.Series([0, 0, 0.646, 0.777, 0.89, 0.512, 0.135, 0.546])

	expected = pandas.Series([0, 0.653])

@patch('widgets._get_git_log')
def test_get_commit_hash(filename_mock):
	test_file = """
	045a5b605b03f566c527f6684586322708525522 045a5b605b03f566c527f6684586322708525522 cdeitrick <cld100@pitt.edu> 1551711670 -0500	checkout: moving from master to version0.2
	045a5b605b03f566c527f6684586322708525522 78db720e4429e60d2821125247c486996d83cc0e Unknown <cld100@pitt.edu> 1551711685 -0500	commit: Update based on pycharm code inspecter
	78db720e4429e60d2821125247c486996d83cc0e d0aa33355336fa3772da8e823660c61296960dfe Unknown <cld100@pitt.edu> 1551713873 -0500	commit: Refactored difference calculation
	d0aa33355336fa3772da8e823660c61296960dfe f086ec9486ea2756f4dd79464c40bfdb02761002 Unknown <cld100@pitt.edu> 1551713984 -0500	commit: Changed Default Clustering Method
	"""
	expected_hash = "f086ec9"

	filename_mock.return_value = test_file
	result_hash = widgets.get_commit_hash()

	assert expected_hash == result_hash

def test_format_linkage_matrix():

	string = """
			left	right	distance	observations
		7	18	0.033999578315858	2
		13	17	0.17508789172405	2
		8	11	0.199140037464566	2
		2	5	0.238709275657774	2
		10	3	0.278982267870099	2
		9	12	0.370131434040108	2
		23	6	0.528725037646305	3
		22	21	0.624301321943297	4
		26	1	0.708258027601832	5
		24	16	0.760211707999897	3
		14	25	0.785856622308224	4
		15	20	0.987877097254146	3
		29	27	1.09399939089444	9
		31	19	1.35849899883873	11
		30	28	1.36249792376227	6
		4	32	1.49918800504887	12
		33	0	2.12548563025534	7
		34	35	4.94288036257192	19
	"""

	expected_table = """
		left	right	distance	observations	clusterId
		7	18	0.033999578315858	2	19
		13	17	0.17508789172405	2	20
		8	11	0.199140037464566	2	21
		2	5	0.238709275657774	2	22
		10	3	0.278982267870099	2	23
		9	12	0.370131434040108	2	24
		23	6	0.528725037646305	3	25
		22	21	0.624301321943297	4	26
		26	1	0.708258027601832	5	27
		24	16	0.760211707999897	3	28
		14	25	0.785856622308224	4	29
		15	20	0.987877097254146	3	30
		29	27	1.09399939089444	9	31
		31	19	1.35849899883873	11	32
		30	28	1.36249792376227	6	33
		4	32	1.49918800504887	12	34
		33	0	2.12548563025534	7	35
		34	35	4.94288036257192	19	36
		"""
	test_table = dataio.import_table(string)
	test_table = widgets.format_linkage_matrix(test_table, 19)
	expected_table = dataio.import_table(expected_table)

	pandas.testing.assert_frame_equal(test_table, expected_table)