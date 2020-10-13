from unittest.mock import patch
from loguru import logger
import pandas
import pytest

from muller import dataio, widgets

@pytest.mark.parametrize(
	"columns, expected",
	[
		(['1', '66', '0', 'X9', 'xc', 'x33', 'col4'], ['1', '66', '0', 'X9', 'x33']),
		( ['Genotype', 0, 1 ,2, 3], [0, 1, 2, 3])
	]
)
def test_get_numeric_columns(columns, expected):
	result = widgets.get_numeric_columns(columns)
	assert result == expected


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


@pytest.mark.parametrize(
	"left,right,index",
	[
		([0, 1, 1, 4, 5], [.23, .14, .13, 0, 0], [0, 1, 2, 3, 4]),
		([0, 1, 0, 0.2, 0], [0, .14, 0, 0, 0], [1, 2, 3]),
		([0, 0, 0, 0, 0], [0, .14, .23, 0, 0], [1, 2]),
		([0, 0, 0, 0, 0], [0, .14, 0, 0, 0], [1]),
		([0, 0, 0, 0, 0], [0, .14, 0, 1, 1], [1, 2, 3, 4]),
	]
)
def test_get_detected_points(left, right, index):
	l = pandas.Series(left)
	r = pandas.Series(right)
	rl, rr = widgets.get_valid_points(l, r, 0.03)
	assert list(rl.index) == list(rr.index)
	assert list(rl.index) == index


def test_get_detected_points_advanced():
	left = pandas.Series([0, 0, 0, 0, 0])
	right = pandas.Series([0, .14, 0, 1, 1])
	result_left, result_right = widgets.get_detected_points(left, right, 0.03, 0.97)
	assert list(result_left.index) == list(result_right.index)
	assert list(result_left.index) == [1]

	left = pandas.Series([0, 0, 0, 0, 0, 1, 1])
	right = pandas.Series([0, 0, 0, .14, .53, 1, 1])
	result_left, result_right = widgets.get_detected_points(left, right, 0.03, 0.97)
	assert list(result_left.index) == list(result_right.index)
	assert list(result_left.index) == [3, 4]
	# Check the `inner` option.
	left = pandas.Series([0, 0, .3, .4, .4, .4, 1, 1])
	right = pandas.Series([0, 0, 0, .1, .1, .1, .2, 1])
	assert [2, 3, 4, 5, 6, 7] == list(widgets.get_detected_points(left, right, .03, inner = False)[0].index)
	assert [2, 3, 4, 5, 6] == list(widgets.get_detected_points(left, right, .03, .97, inner = False)[0].index)
	assert [3, 4, 5, 6, 7] == list(widgets.get_detected_points(left, right, .03, inner = True)[0].index)
	assert [3, 4, 5] == list(widgets.get_detected_points(left, right, .03, .97, inner = True)[0].index)

def test_get_detected_points_inner():
	left = pandas.Series([0, 0, 0, 0,   0,    0, 0.085, 0.001, 0.005])
	right = pandas.Series([0,0, 0,   0,   0,  0,0.05, 0.55, 0.5 ])
	l,r = widgets.get_valid_points(left, right, dlimit = 0.03, inner = True)

	assert l.tolist() == [0.085]
	assert r.tolist() == [0.05]


def test_get_valid_points_simple():
	left = pandas.Series([0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 0])
	right = pandas.Series([0, 0, 0, .1, .2, .3, .3, .3, .3, 0, 0, 0])

	result_left, result_right = widgets.get_valid_points(left, right, 0.03)
	assert result_left.tolist() == [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]
	assert result_right.tolist() == [0, 0, .1, .2, .3, .3, .3, .3, 0, 0]

	result_left, result_right = widgets.get_valid_points(left, right, 0.03, 0.97)
	assert result_left.tolist() == [.1, .2, .3, .4, .5, .6, .7, .8, .9]
	assert result_right.tolist() == [0, 0, .1, .2, .3, .3, .3, .3, 0]

	expected = pandas.DataFrame({
		'left':  [.3, .4, .5, .6, .7, .8],
		'right': [.1, .2, .3, .3, .3, .3],
	}, index = range(3, 9))
	result_left, result_right = widgets.get_valid_points(left, right, 0.03, 0.97, inner = True)
	assert result_left.tolist() == [.3, .4, .5, .6, .7, .8]
	assert result_right.tolist() == [.1, .2, .3, .3, .3, .3]


def test_get_valid_points_complex():
	left = pandas.Series([0.00, 0.00, 0.000, 0.00, 0.00, 0.263, 0.07, 0.081, 0.069, 0.042])
	right = pandas.Series([0.00, 0.00, 0.170, 0.55, 0.947, 1.00, 1.00, 1.00, 1.00, 1.00])

	expected_left = [0.000, 0.00, 0.00, 0.263, 0.07, 0.081, 0.069, 0.042]
	expected_right = [0.170, 0.55, 0.947, 1.00, 1.00, 1.00, 1.00, 1.00]
	result_left, result_right = widgets.get_valid_points(left, right, dlimit = 0.03)
	assert result_left.tolist() == expected_left
	assert result_right.tolist() == expected_right

	switched_result_left, switched_result_right = widgets.get_valid_points(right, left, 0.03)
	assert switched_result_left.tolist() == expected_right
	assert switched_result_right.tolist() == expected_left

	expected_left = [0.263, 0.07, 0.081, 0.069, 0.042]
	expected_right = [1.00, 1.00, 1.00, 1.00, 1.00]
	result_left, result_right = widgets.get_valid_points(left, right, 0.03, inner = True)
	assert result_left.tolist() == expected_left
	assert result_right.tolist() == expected_right

	result_left, result_right = widgets.get_valid_points(left, right, 0.03, 0.97, inner = True)
	assert result_left.tolist() == [] and result_right.tolist() == []


@patch('muller.widgets._get_git_log')
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





@pytest.mark.parametrize(
	"series,expected",
	[
		([0, 0, 0, 1, 1, 1], True),
		([0, 1, 0, 1, 0, 1], True),
		([0, .2, 1, 1, 1], False)
	]
)
def test_fixes_imediately(series, expected):
	s = pandas.Series(series)
	assert widgets.fixed_immediately(s, 0.03, 0.97) == expected


@pytest.mark.parametrize(
	"series,expected",
	[
		([0, 0, 0, 1, 1, 1], True),
		([0, 1, 0, 1, 0, 1], True),
		([0, .2, 1, 1, 1], True),
		([0, .1, .2, .3, .4, .5], False)
	]
)
def test_fixed(series, expected):
	s = pandas.Series(series)
	assert widgets.fixed(s, 0.97) == expected

@pytest.mark.parametrize(
	"series, expected",
	[
		([0.1, 0.4, 0.3, 1, 0.97 , 1], (3,5)),
		([0.2, 1, 0.2, 0.98, 0.1], (1,3)),
		([0.1, 0.2, 0.3, 0.4, 0.5], None),
		([0.1, 0.4, 0.3, 1, 0.5, 0.2], (3,3))
	]
)
def test_find_boundaries_fixed(series, expected):
	s = pandas.Series(series)
	result = widgets.find_boundaries_fixed(s, 0.97)

	assert result == expected


@pytest.mark.parametrize(
	"left, right, expected",
	[
		([0, 0, 0.261, 1.000, 1.000, 1.000, 1.000],[0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.910], [0,0,0,0,0,1,1]),
		([0, 0, 0.261, 1.000, 0.000, 0.200, 0.200],[0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.910], [0,0,0,0,0,0,0]),
	]
)
def test_get_overlap_regions(left, right, expected):
	result = widgets.get_overlap_regions(left, right, 0.9)
	# Convert to bool for safety.
	assert result.tolist() == [bool(i) for i in expected]

@pytest.mark.parametrize(
	"values, expected",
	[
		(pandas.Series([4,7,9,11]), [4,7,9,11]),
		([1,88,4,88], [1,88,4,88]),
		('string1', ['s', 't', 'r', 'i', 'n', 'g', '1'])
	]
)
def test_coerce_to_list(values, expected):
	result = widgets._coerce_to_list(values)
	assert result == expected

@pytest.mark.parametrize(
	"values, expected",
	[
		([0.000, 0.000, 0.261, 1.000, 1.000, 1.000, 1.000], 3),
		([0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.910], 5),
		([.1, .1, .1, .1], None),
		([1, .1, .1, .1], 0)
	]

)
def test_get_first_fixed_timepoint(values, expected):
	result = widgets.get_first_fixed_timepoint(values, 0.9)
	assert result == expected

@pytest.mark.parametrize(
	"values, expected",
	[
		([0.000, 0.000, 0.261, 1.000, 1.000, 1.000, 1.000], pandas.Series([0.261], index = [2])),
		([0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.910], pandas.Series([0.525, 0.454], index = [3, 4])),
		([0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.810], pandas.Series([0.525, 0.454, 0.810], index = [3,4,6])),
		([0.000, 0.000, 1.000, 1.000, 1.000, 1.000, 1.000], pandas.Series([]))
	]
)
def test_get_intermediate(values, expected):
	result = widgets.get_intermediate(values, 0.03, 0.9)

	# Since pandas likes to return a series of bool values when comparing items rather than a scalar result,
	# Let's check the values and index directly.
	assert result.tolist() == expected.tolist()
	assert list(result.index) == list(expected.index)

@pytest.mark.parametrize(
	"values, expected",
	[
		([0.000, 0.000, 0.261, 1.000, 1.000, 1.000, 1.000], pandas.Series([1,1,1,1], index = [3,4,5,6])),
		([0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.910], pandas.Series([0.911, 0.910], index = [5,6])),
		([0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.810], pandas.Series([0.911], index = [5])),
		([0.000, 0.000, 0.860, 0.000, 0.000, 0.000, 0.000], pandas.Series([]))
	]
)
def test_get_fixed(values, expected):
	result = widgets.get_fixed(values, 0.9)

	# Since pandas likes to return a series of bool values when comparing items rather than a scalar result,
	# Let's check the values and index directly.
	assert result.tolist() == expected.tolist()
	assert list(result.index) == list(expected.index)

@pytest.mark.parametrize(
	"values, expected",
	[
		([0.000, 0.000, 0.261, 1.000, 1.000, 1.000, 1.000], pandas.Series([0,0], index = [0,1])),
		([0.000, 0.000, 0.000, 0.525, 0.454, 0.911, 0.910], pandas.Series([0, 0,0], index = [0,1,2])),
		([0.000, 0.000, 0.000, 0.525, 0.020, 0.911, 0.810], pandas.Series([0,0,0,0.020], index = [0,1,2,4])),
		([1.000, 1.000, 0.860, 1.000, 1.000, 1.000, 1.000], pandas.Series([]))
	]
)
def test_get_undetected(values, expected):
	result = widgets.get_undetected(values, 0.03)

	# Since pandas likes to return a series of bool values when comparing items rather than a scalar result,
	# Let's check the values and index directly.
	assert result.tolist() == expected.tolist()
	assert list(result.index) == list(expected.index)

@pytest.mark.parametrize(
	"elements, size, expected",
	[
		(3, 3, 1),
		(4, 2, 6),
		(6, 3, 20)
	]
)
def test_calculate_total_number_of_combinations(elements, size, expected):

	result = widgets.calculate_number_of_combinations(elements, size)
	assert result == expected


