import pytest
from muller.inheritance import scoring
import pandas
from loguru import logger
@pytest.mark.parametrize(
	"left,right,expected",
	[	# nested, unnested, expected_score
		([0, 0.2, 0.3, 0.4, 0.5], [0, 0.1, 0.2, 0.3, 0.4], 1),
		([0.3, 0.3, 0.3, 0.3, 0.3], [0, 0.1, 0.2, 0.3, 0.4], 0),
		([0.3, 0.3, 0.3, 0.3, 0.3], [0.3, 0.4, 0.5, 0.6, 0.7], -1)
	]
)
def test_summation_score(left, right, expected):
	left_series = pandas.Series(left)
	right_series = pandas.Series(right)

	result  = scoring.calculate_subtractive_score(left_series, right_series, 0.1)

	assert result == expected


@pytest.mark.parametrize(
	"left,right,expected",
	[	# nested, unnested, expected_score
		([0, 0.2, 0.3, 0.4, 0.5], [0, 0.1, 0.2, 0.3, 0.4], 2),
		([0.3, 0.3, 0.3, 0.3, 0.3], [0, 0.1, 0.2, 0.1, 0.2], 0),
		([0.3, 0.2, 0.1, 0.0, 0.0], [0.3, 0.4, 0.5, 0.6, 0.7], -2)
	]
)
def test_calculate_derivative_score(left, right, expected):
	left_series = pandas.Series(left)
	right_series = pandas.Series(right)

	result = scoring.calculate_derivative_score(left_series, right_series, 0.05, 0.01)
	assert result == expected


def test_calculate_subtractive_score():
	pass

