from typing import List

import pandas
import pytest

from inheritance import checks


@pytest.mark.parametrize("left,right,expected",
	[
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], False),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .93, 0], False),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], True),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, 1.0, 0], True),
		([0,0,0,0.403,0.489,0.057,0.08], [0,0,0,0,0,0.2675,0.326], False)
	])
def test_additive_check(left: List[float], right: List[float], expected: bool):
	# check if two series consistently sum to greater than 1.
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = checks.check_additive_background(left, right, 1.03, 1.15)

	assert expected == result


@pytest.mark.parametrize("left,right,expected",
	[
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], False),
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .4, .3], False),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .93, 0], True),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], True),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, 1.0, 0], True)
	])
def test_subtractive_check(left: List[float], right: List[float], expected: bool):
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = checks.check_subtractive_background(left, right, -.03, -0.15)

	assert expected == result

@pytest.mark.parametrize("left,right,expected",
	[
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .3, .3], True),
		([0, .1, .1, .2, .2, .3, .3], [0, .1, .1, .2, .2, .4, .3], True),
		([1, .9, .8, .7, .6, .5, .4], [0, .1, .2, .3, .4, .5, .6], False),
		([0, .1, .1, .1, .1, .1, .1], [0, 0, 0, 0, 0, .94, .94], True),
		([0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1.0, 0], False)
	])
def test_check_derivative_background(left: List[float], right: List[float], expected: bool):
	left = pandas.Series(left)
	right = pandas.Series(right)

	result = checks.check_derivative_background(left, right, -.03)

	assert expected == (result > 0)