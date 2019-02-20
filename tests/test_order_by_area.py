import pytest
import pandas
from inheritance import order_by_area
from typing import Tuple
@pytest.fixture
def full_overlap()->Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, 0, 0, 0, 1, 1])
	right=  pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right

@pytest.fixture
def no_overlap()->Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0,.3, .4, 0, 0, 0])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right

@pytest.fixture
def partial_overlap()->Tuple[pandas.Series, pandas.Series]:
	left = pandas.Series([0, 0, .1, .2, 1, 1])
	right = pandas.Series([0, 0, 0, 0, .2, .3])
	return left, right


def test_calculate_overlapping_area(full_overlap, no_overlap, partial_overlap):
	left, right = full_overlap
	right_area = order_by_area.area_of_series(right)
	overlap_area = order_by_area.calculate_common_area(left, right, 0.03)
	assert overlap_area == right_area

	left, right = no_overlap
	overlap_area = order_by_area.calculate_common_area(left, right, 0.03)
	assert overlap_area == 0

	left, right = partial_overlap
	overlap_area = order_by_area.calculate_common_area(left, right, 0.03)
	assert overlap_area == right_area

	left = pandas.Series([0.01,0.279,0.341,0.568,0.708,0.913,0.756,0.455,0.399,0.13,0.041])
	right = pandas.Series([0,0,0,0,0,0.247,0.388,0.215,0.403,0.141,0.028])
	overlap_area = order_by_area.calculate_common_area(left, right, 0)
	assert overlap_area == order_by_area.area_of_series(right)
def test_calculate_noncommon_area(full_overlap, partial_overlap, no_overlap):
	left, right = full_overlap
	expected_area = abs(order_by_area.area_of_series(right) - order_by_area.area_of_series(left))
	nonoverlap_area = order_by_area.calculate_noncommon_area(left, right, 0.03)
	assert nonoverlap_area == expected_area

	left, right = partial_overlap
	expected_area = order_by_area.area_of_series(pandas.Series([.1, .2, .8, .7]))
	nonoverlap_area = order_by_area.calculate_noncommon_area(left, right, 0.03)
	assert nonoverlap_area == expected_area

	left, right = no_overlap
	expected_area = order_by_area.area_of_series(left) + order_by_area.area_of_series(right)
	noncommon_area = order_by_area.calculate_noncommon_area(left, right, .03)
	assert expected_area == noncommon_area