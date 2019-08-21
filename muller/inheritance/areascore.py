import math
from typing import Set, Union

import pandas


def area_of_series(series: pandas.Series) -> float:
	""" Calculates the area of a discrete series."""
	total = series.sum()
	return total


def calculate_common_area(left: pandas.Series, right: pandas.Series) -> float:
	""" Calculates |X and Y|"""
	overlap_index = get_overlap(left, right)
	if not overlap_index:
		overlap = 0
	else:
		overlap = area_of_series(left.where(left < right, right))
	return overlap


def calculate_area_difference(left: pandas.Series, right: pandas.Series):
	left_area = area_of_series(left)
	right_area = area_of_series(right)
	return left_area - right_area


# noinspection PyTypeChecker
def get_overlap(left: pandas.Series, right: pandas.Series) -> Set[int]:
	detected_left = left
	detected_right = right

	overlap_index = set(detected_left.index) & set(detected_right.index)
	return overlap_index


def X_or_Y(left: Union[float, pandas.Series], right: Union[float, pandas.Series], common: float = None) -> float:
	""" Calculated the XOR area between `left` and `right`
		Parameters
		----------
		left, right: pandas.Series
			The two series to calculate the exclusive area on.
		left, right,common: float
			Precomputed values for the area of `left`, `right`, and x_and_y.
	"""
	if common is None:
		area_left = area_of_series(left)
		area_right = area_of_series(right)
		common_area = calculate_common_area(left, right)
	else:
		area_left, area_right, common_area = left, right, common

	return area_left + area_right - common_area


def jaccard_distance(left: Union[float, pandas.Series], right: Union[float, pandas.Series]) -> float:
	""" Calulated the jaccard distance between two series.
		Parameters
		----------
		left, right: pandas.Series
			The two genotypes to calculate the jaccard distance on.
		left, right: float
			precomputed x_or_y and x_and_y, respectively.
	"""
	if isinstance(left, float):
		x_or_y, x_and_y = left, right
	else:
		x_or_y = X_or_Y(left, right)
		x_and_y = X_and_Y(left, right)
	return (x_or_y - x_and_y) / x_or_y


def jaccard_subset(left: Union[float, pandas.Series], right: Union[float, pandas.Series]) -> float:
	"""
		Calculates the jaccard distance assuming `right` is a subset of `left`.
	Parameters
	----------
	left, right: pandas.Series
	left, right: float
	"""
	if isinstance(left, float):
		left_area, right_area = left, right
	else:
		left_area = area_of_series(left)
		right_area = area_of_series(right)
	return (left_area - right_area) / left_area


def is_subset_legacy(left: pandas.Series, right: pandas.Series) -> bool:
	""" Tests if `right` is a subset of `left`. Deprecated since it needs to be faster."""
	actual_jaccard_distance = jaccard_distance(left, right)
	subset_jaccard_distance = jaccard_subset(left, right)
	# TODO: incorporate variance into this
	return math.isclose(actual_jaccard_distance, subset_jaccard_distance, abs_tol = 0.1)


def is_subset(left: Union[float, pandas.Series], right: Union[float, pandas.Series], common: float = None):
	"""
		Estimates whether `right` is a subset of `left`. This implementation should be a little faster than the legacy method.
	Parameters
	----------
	left, right: pandas.Series
	left,right,common: float
		precomuted areas
	"""
	if common is None:
		left_area = area_of_series(left)
		right_area = area_of_series(right)
		x_and_y = calculate_common_area(left, right)
	else:
		left_area, right_area, x_and_y = left, right, common
	x_or_y = X_or_Y(left_area, right_area, x_and_y)

	actual_jaccard_distance = jaccard_distance(x_or_y, x_and_y)
	subset_jaccard_distance = jaccard_subset(left_area, right_area)  # Jaccard distance assuming `left` is a subset of `right`
	return math.isclose(actual_jaccard_distance, subset_jaccard_distance, abs_tol = 0.1)


X_and_Y = calculate_common_area
