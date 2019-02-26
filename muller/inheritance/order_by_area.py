import math
from typing import Set

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


def X_or_Y(left: pandas.Series, right: pandas.Series) -> float:
	area_left = area_of_series(left)
	area_right = area_of_series(right)
	common_area = calculate_common_area(left, right)

	return area_left + area_right - common_area


def jaccard_distance(left: pandas.Series, right: pandas.Series) -> float:
	x_or_y = X_or_Y(left, right)
	return (x_or_y - X_and_Y(left, right)) / x_or_y


def jaccard_subset(left: pandas.Series, right: pandas.Series) -> float:
	left_area = area_of_series(left)
	right_area = area_of_series(right)
	return (left_area - right_area) / left_area


def is_subset(left: pandas.Series, right: pandas.Series) -> bool:
	""" Tests if `right` is a subset of `left`"""
	actual_jaccard_distance = jaccard_distance(left, right)
	subset_jaccard_distance = jaccard_subset(left, right)

	# TODO: incorporate variance into this
	return math.isclose(actual_jaccard_distance, subset_jaccard_distance, abs_tol = 0.1)


X_and_Y = calculate_common_area

if __name__ == "__main__":
	left = pandas.Series([0, .1, .2, .3, .4])
	right = pandas.Series([0, 0, .1, .2, .3])

	print(is_subset(left, right))
	print(is_subset(right, left))
