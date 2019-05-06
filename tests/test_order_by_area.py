from typing import Tuple

import pandas
import pytest

from muller.inheritance import order_by_area


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


def test_calculate_overlapping_area(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	right_area = order_by_area.area_of_series(right)
	overlap_area = order_by_area.calculate_common_area(left, right)
	assert overlap_area == right_area

	left, right = no_overlap
	overlap_area = order_by_area.calculate_common_area(left, right)
	assert overlap_area == 0

	left, right = partial_overlap
	overlap_area = order_by_area.calculate_common_area(left, right)
	assert overlap_area == right_area

	left = pandas.Series([0.01, 0.279, 0.341, 0.568, 0.708, 0.913, 0.756, 0.455, 0.399, 0.13, 0.041])
	right = pandas.Series([0, 0, 0, 0, 0, 0.247, 0.388, 0.215, 0.399, 0.13, 0.028])
	overlap_area = order_by_area.calculate_common_area(left, right)
	assert overlap_area == order_by_area.area_of_series(right)


def test_is_subset(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	assert order_by_area.is_subset(left, right)
	assert False == order_by_area.is_subset(right, left)

	left, right = partial_overlap
	assert True == order_by_area.is_subset(left, right)
	assert False == order_by_area.is_subset(right, left)
	left, right = no_overlap
	assert False == order_by_area.is_subset(left, right)
	assert False == order_by_area.is_subset(right, left)
