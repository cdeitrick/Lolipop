from typing import List

import pandas
import pytest

from muller.inheritance import scoring


@pytest.mark.parametrize("left,right,expected",
	[
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], 0),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .93, 0], -1),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], -1),
		([0, .1, .1, .1, .1, .16, .1], [0, 0, 0, 0, 0, 1.0, 0], -1),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], [0, 0, 0, 0, 0, 0.2675, 0.326], 0),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], [0, 0, 0, 0.2, 0.2, 0, 0], 1)
	])
def test_calculate_additive_score(left: List[float], right: List[float], expected: int):
	# check if two series consistently sum to greater than 1.
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = scoring.calculate_greater_score(left, right, 0.03)

	assert result == expected


@pytest.mark.parametrize("left,right,expected",
	[
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], 0),
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .4, .3], 0),
		([1, .9, .8, .7, .6, .5, .4], [0, .1, .2, .3, .4, .5, .6], -2),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], 0),
		([0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1.0, 0], 0)
	])
def test_calculate_derivative(left: List[float], right: List[float], expected: int):
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = scoring.calculate_derivative_score(left, right, 0.03, 0.03)

	assert result == expected
