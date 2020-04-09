from pathlib import Path
from typing import *
import pytest
from muller.clustering.metrics import distance_calculator


@pytest.mark.parametrize(
	"left, right, expected",
	[
		([0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00], [0.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00], "onlyFixed"),
		([0.00, 0.00, 0.00, 1.00, 1.00, 1.00, 1.00], [0.00, 0.00, 0.00, 0.52, 0.45, 0.91, 0.91], "partiallyFixed"),
		([0.00, 0.01, 0.26, 1.00, 1.00, 1.00, 1.00], [0.00, 0.00, 0.00, 0.52, 0.45, 0.91, 0.91], "bothFixed"),
		([0.00, 0.01, 0.26, 1.00, 1.00, 1.00, 1.00], [0.00, 0.00, 0.00, 0.18, 0.17, 0.23, 0.24], "oneFixed"),
		([0.00, 0.00, 0.00, 0.11, 0.00, 0.11, 0.12], [0.00, 0.00, 0.00, 0.18, 0.17, 0.23, 0.24], "notFixed"),
	]
)
def test_categorize_series(left, right, expected):
	result = distance_calculator.get_pair_category(left, right, dlimit = 0.03, flimit = 0.90)
	assert result == expected