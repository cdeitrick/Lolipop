import pandas
import pytest

from muller.clustering.metrics import distance_methods


@pytest.mark.parametrize("left,right,expected",
	[
		([1, 7, 6, 2, 4], [8, 6, 3, 5, 6], 8.48528),
		([0, 0.0, 0.0, 0.273, 0.781, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.345, 0.833, 0.793], 0.579105),
		([0, 0.0, 0.261, 1.0, 1.0, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0], 1.4381)
	]
)
def test_minkowski_distance(left, right, expected: float):
	left = pandas.Series(left)
	right = pandas.Series(right)
	result = distance_methods.minkowski_distance(left, right, 2)
	assert pytest.approx(result, rel = 1E-4) == expected


@pytest.mark.parametrize("left,right,expected",
	[
		([0, 0.0, 0.0, 0.273, 0.781, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.345, 0.833, 0.793], 1 - 0.945696),
		([1, 2, 2, 1, 5, 1, 2, 6, 2], [12, 4, 1, 4, 5, 1, 4, 1, 3], 1 + 0.274663),
		([0, 0.0, 0.261, 1.0, 1.0, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0], 1 - 0.539221)
	]
)
def test_pearson_distance(left, right, expected: float):
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = distance_methods.pearson_correlation_distance(left, right, adjusted = False)
	assert pytest.approx(result, rel = 1E-4) == expected


@pytest.mark.parametrize("left,right",
	[
		([0, 0.0, 0.0, 0.273, 0.781, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.345, 0.833, 0.793]),
		([.1, .2, .2, .1, .5, .1, .2, .6, .2], [.12, .4, .1, .4, .5, .1, .4, .1, .3]),
		([0, 0.0, 0.261, 1.0, 1.0, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0])
	]
)
def test_calculate_distance(left, right):
	left = pandas.Series(left)
	right = pandas.Series(right)
	pearson = distance_methods.pearson_correlation_distance(left, right)
	minkowski = distance_methods.minkowski_distance(left, right)
	jaccard = distance_methods.jaccard_distance(left, right)
	binomial = distance_methods.binomial_distance(left, right)

	assert pearson == distance_methods.calculate_distance(left, right, 'pearson')
	assert minkowski == distance_methods.calculate_distance(left, right, 'minkowski')
	assert jaccard == distance_methods.calculate_distance(left, right, 'jaccard')
	assert binomial == distance_methods.calculate_distance(left, right, 'binomial')
