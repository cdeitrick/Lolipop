from typing import Tuple

import pandas
import pytest
from loguru import logger

from muller.inheritance import areascore
from ..filenames import real_tables


#################################################################################
############################### Fixtures ########################################
#################################################################################

@pytest.fixture
def full_overlap() -> Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, 0, 0, 0, 1, 1])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


@pytest.fixture
def no_overlap() -> Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, .3, .4, 0, 0, 0])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


@pytest.fixture
def partial_overlap() -> Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, 0, .1, .2, 1, 1])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


#################################################################################
############################### Small Tests #####################################
#################################################################################

def test_calculate_overlapping_area(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	right_area = areascore.area_of_series(right)
	overlap_area = areascore.calculate_common_area(left, right)
	assert pytest.approx(overlap_area, abs = 0.001) == right_area

	left, right = no_overlap
	overlap_area = areascore.calculate_common_area(left, right)
	assert pytest.approx(overlap_area, abs = 0.001) == 0

	left, right = partial_overlap
	overlap_area = areascore.calculate_common_area(left, right)
	assert pytest.approx(overlap_area, abs = 0.001) == right_area

	left = pandas.Series([0.01, 0.279, 0.341, 0.568, 0.708, 0.913, 0.756, 0.455, 0.399, 0.13, 0.041])
	right = pandas.Series([0, 0, 0, 0, 0, 0.247, 0.388, 0.215, 0.399, 0.13, 0.028])
	overlap_area = areascore.calculate_common_area(left, right)
	assert pytest.approx(overlap_area) == areascore.area_of_series(right)


@pytest.mark.parametrize(
	"left, right, expectedforward, expectedreverse",
	[  # Full overlap
		([0, 0, 0, 0, 1, 1], [0, 0, 0, 0, .2, .3], True, False),
		# Partial overlap
		([0, 0, .1, .2, 1, 1], [0, 0, 0, 0, .2, .3], True, False),
		# No overlap
		([0, .3, .4, 0, 0, 0], [0, 0, 0, 0, .2, .3], False, False),
		# Other
		([0, 1, 0, 1, 0], [0, 0,.5, 0, 0], False, False)
	]
)
def test_is_subset(left, right, expectedforward, expectedreverse):
	""" `expectedforward` refers to whether right is a subset of left,
		`expectedreverse` refers to whether left is a subset of right.
	"""
	left = pandas.Series(left)
	right = pandas.Series(right)
	assert areascore.is_subset(left, right) == expectedforward
	assert areascore.is_subset(right, left) == expectedreverse


"""
def test_is_subset(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	assert areascore.is_subset(left, right)
	assert False == areascore.is_subset(right, left)

	left, right = partial_overlap
	assert True == areascore.is_subset(left, right)
	assert False == areascore.is_subset(right, left)
	left, right = no_overlap
	assert False == areascore.is_subset(left, right)
	assert False == areascore.is_subset(right, left)
"""


@pytest.mark.parametrize(
	"series, expected",
	[
		([0, 0.2, 0.3, 0.4, 0.5], 1.4),
		([0.3, 0.3, 0.3, 0.3, 0.3], 1.5),
		([0, .1, .1, .2, .2, .3, .3], 1.2),
		([0, 0, 0, 0.403, 0.489, 0.057, 0.08], 1.029),
		([0, 0, 0, 0, 0], 0),
		([0, 1, 0, .1, 0], 0.55)
	]
)
def test_area_of_series(series, expected):
	s = pandas.Series(series)
	result = areascore.area_of_series(s)

	assert pytest.approx(result, expected)


def test_area_of_series_real():
	filename = list(real_tables.values())[0]
	table_genotype = pandas.read_excel(filename, sheet_name = "genotype").set_index('Genotype')
	for index, row in table_genotype.iterrows():
		logger.debug(index)
		expected_area = areascore.calculate_area(row)
		actual_area = areascore.area_of_series(row)
		assert pytest.approx(actual_area, abs = .001) == expected_area
