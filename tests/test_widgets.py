from unittest.mock import patch

import pandas

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


def test_get_valid_points():
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


def test_get_commit_hash():
	test_file = """
	045a5b605b03f566c527f6684586322708525522 045a5b605b03f566c527f6684586322708525522 cdeitrick <cld100@pitt.edu> 1551711670 -0500	checkout: moving from master to version0.2
	045a5b605b03f566c527f6684586322708525522 78db720e4429e60d2821125247c486996d83cc0e Unknown <cld100@pitt.edu> 1551711685 -0500	commit: Update based on pycharm code inspecter
	78db720e4429e60d2821125247c486996d83cc0e d0aa33355336fa3772da8e823660c61296960dfe Unknown <cld100@pitt.edu> 1551713873 -0500	commit: Refactored difference calculation
	d0aa33355336fa3772da8e823660c61296960dfe f086ec9486ea2756f4dd79464c40bfdb02761002 Unknown <cld100@pitt.edu> 1551713984 -0500	commit: Changed Default Clustering Method
	"""
	expected_hash = "f086ec9"

	with patch('muller.widgets._get_git_log') as filename_mock:
		filename_mock.return_value = test_file
		result_hash = widgets.get_commit_hash()

	assert expected_hash == result_hash
