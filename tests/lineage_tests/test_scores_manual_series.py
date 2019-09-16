from typing import List

import pandas
import pytest

from muller.inheritance import scoring


@pytest.fixture
def scorer() -> scoring.Score:
	return scoring.Score(0.05, 0.03, 0.97)


@pytest.fixture
def legacy_scorer() -> scoring.LegacyScore:
	return scoring.LegacyScore(0.05, 0.03, 0.97)


@pytest.mark.parametrize(
	"left,right,expected",
	[  # nested, unnested, expected_score
		([0, 0.2, 0.3, 0.4, 0.5], [0, 0.1, 0.2, 0.3, 0.4], 1),
		([0.3, 0.3, 0.3, 0.3, 0.3], [0, 0.1, 0.2, 0.3, 0.4], 0),
		([0.3, 0.3, 0.3, 0.3, 0.3], [0.3, 0.4, 0.5, 0.6, 0.7], -1),
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], 0),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .93, 0], 0),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], 0),
		([0, .1, .1, .1, .1, .16, .1], [0, 0, 0, 0, 0, 1.0, 0], 0),
		# ([0, .1, .2, .3, .4, .16, .1], [0, 0, 0, 0, 0, 1.0, 0], -1),
		# ([0, .2, .2, .2, .4, .2, .2], [0, 0, 0, 0, 0, 1.0, 0], -1),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], [0, 0, 0, 0, 0, 0.2675, 0.326], 0),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], [0, 0, 0, 0.2, 0.2, 0, 0], 1)
	]
)
def test_greater_score(scorer, left, right, expected):
	left_series = pandas.Series(left)
	right_series = pandas.Series(right)

	result = scorer.calculate_score_greater(left_series, right_series)

	assert result == expected


@pytest.mark.parametrize("left,right,expected",
	[
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], 0),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], [0, 0, 0, 0, 0, 0.2675, 0.326], 0),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], 0),
		([0, .1, .1, .1, .1, .1, .1], [0, .9, .9, 1.0, .5, .93, .9], 0),
		([0, .1, .2, .3, .3, .2, .1], [0, .9, .9, 1.0, .5, .93, .9], 1),
		([0, .1, .1, .3, .5, .5, .2], [0.2, 0.9, 0.85, 0.9, .95, 1.0, 0.9], 1),
		([0, 0.5, 0.5, 0.403, 0.489, 0.05, 0.05], [0, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7], 1)
	])
def test_calculate_above_fixed_score(scorer, left: List[float], right: List[float], expected: int):
	# check if two series consistently sum to greater than 1.
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = scorer.calculate_score_above_fixed(left, right)

	assert result == expected


@pytest.mark.parametrize(
	"left,right,expected",
	[  # nested, unnested, expected_score
		([0, 0.2, 0.3, 0.4, 0.5], [0, 0.1, 0.2, 0.3, 0.4], 2),
		([0.3, 0.3, 0.3, 0.3, 0.3], [0, 0.1, 0.2, 0.1, 0.2], 0),
		([0.3, 0.2, 0.1, 0.0, 0.0], [0.3, 0.4, 0.5, 0.6, 0.7], -2),
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], 0),
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .4, .3], 0),
		([1, .9, .8, .7, .6, .5, .4], [0, .1, .2, .3, .4, .5, .6], -2),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], 0),
		([0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1.0, 0], 0)
	]
)
def test_calculate_derivative_score(legacy_scorer, left, right, expected):
	left_series = pandas.Series(left)
	right_series = pandas.Series(right)

	result = legacy_scorer.calculate_derivative_score(left_series, right_series)
	assert result == expected


@pytest.mark.parametrize(
	"left, right, expected",
	[
		([0.1, 0.1, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.2, 0.0, 0.0], -2),
		([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],0),
		([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2],2),
		([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [0.0, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],2),
		([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0, 0.0],2),
	]
)
def test_calculate_area_score(scorer, left, right, expected):
	l = pandas.Series(left)
	r = pandas.Series(right)
	result = scorer.calculate_score_area(l, r)

	assert result == expected
