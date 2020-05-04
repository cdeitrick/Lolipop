from pathlib import Path

import pandas
import pytest
from loguru import logger
from muller import dataio
from muller.clustering.metrics import DistanceCalculator, distance_calculator, distance_methods
import math
from tests import filenames


@pytest.fixture
def b1_data() -> pandas.DataFrame:
	f = filenames.real_tables["B1"]
	t = dataio.import_table(f, sheet_name = 'trajectory', index = 'Trajectory')
	t.index = [str(i) for i in t.index]
	return t


@pytest.fixture
def pairwise_values(b1_data):
	dc = DistanceCalculator(0.03, 0.97, 'binomial')
	pairwise_values = dc.run(b1_data)
	return pairwise_values


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
		([0, 0.0, 0.0, 0.273, 0.781, 1.0, 1.0], [0, 0.0, 0.0, 0.0, 0.345, 0.833, 0.793],
		1 - 0.945696),
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


@pytest.mark.parametrize("name,closest",
	[
		("4", "8"),
		("2", "3"),
		("17", "9"),
		("5", "20")
	]
)
def test_binomial_distance(pairwise_values, name: str, closest: str):
	logger.info(pairwise_values.keys())
	candidates = (i for i in pairwise_values.items() if (name in i[0] or int(name) in i[0]))
	pair, value = min(candidates, key = lambda s: s[1])
	assert name in pair and closest in pair


@pytest.mark.parametrize("left,right,expected",
	[
		([0, 0, 0, 1, 1, 1, 1, 0], [0, 0, 0, 1, 1, 1, 1, 0], 0),
		([0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0], 0),
		([0, 0, 0, 1, 1, 1, 1, 0], [0, 0, 0, 1, 1, 1, 1, 1], 1),
		([0, 0, 0, 1, 1, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1], 1),
		([0, 0, 0, 1, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1], 1),
		([0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 1, 1], 1),
		([0, 0, 0, 1, 1, 1, 1], [0, 0, 0, 1, 1, 1, 1], 0)
	]
)
def test_fixed_overlap(left, right, expected):
	left = pandas.Series(left)
	right = pandas.Series(right)
	result = distance_calculator.fixed_overlap(left, right, 0.97)
	if math.isnan(result):
		result = 1
	assert expected == result
