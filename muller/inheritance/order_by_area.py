from typing import Set

import pandas


def area_of_series(series: pandas.Series) -> float:
	""" Calculates the area of a discrete series."""
	total = series.sum()
	return total


def calculate_noncommon_area(left: pandas.Series, right: pandas.Series) -> float:
	""" Calculates the area not contained by bboth `left` and `right`"""
	difference = abs(left - right)
	noncommon_area = area_of_series(difference)

	return noncommon_area


def calculate_non_overlapping_area(left: pandas.Series, right: pandas.Series, dlimit) -> float:
	# Area of the timepoionts that do no overlap
	overlap_index = get_overlap(left, right, dlimit)
	nonoverlap_left = left.loc[~left.index.isin(overlap_index)]
	nonoverlap_right = right.loc[~right.index.isin(overlap_index)]
	area_left = area_of_series(nonoverlap_left)
	area_right = area_of_series(nonoverlap_right)
	return area_left + area_right


def calculate_common_area(left: pandas.Series, right: pandas.Series, dlimit: float) -> float:
	# left = left.reset_index(drop = True)
	# right = right.reset_index(drop = True)
	overlap_index = get_overlap(left, right, dlimit)
	if not overlap_index:
		overlap = 0
	else:
		overlap = area_of_series(left.where(left < right, right))
	return overlap

def calculate_area_difference(left: pandas.Series, right:pandas.Series):
	left_area = area_of_series(left)
	right_area = area_of_series(right)
	return left_area - right_area
# noinspection PyTypeChecker
def get_overlap(left: pandas.Series, right: pandas.Series, dlimit: float) -> Set[int]:
	detected_left = left[left > dlimit]
	detected_right = right[right > dlimit]

	overlap_index = set(detected_left.index) & set(detected_right.index)
	return overlap_index


# noinspection PyTypeChecker
def compare_to_table(left: pandas.Series, table: pandas.DataFrame, dlimit: float) -> pandas.DataFrame:
	area = area_of_series(left)
	axis = 1
	common_area = table.apply(calculate_common_area, axis = axis, args = (left, dlimit)) / area
	noncommon_area = table.apply(calculate_noncommon_area, axis = axis, args = (left, dlimit)) / area
	nonoverlapping_area = table.apply(calculate_non_overlapping_area, axis = axis, args = (left, dlimit)) / area

	df = pandas.concat([common_area, noncommon_area, nonoverlapping_area], axis = 1)
	df.columns = ['commonArea', 'noncommonArea', 'nonoverlappingArea']
	df = df.sort_values(by = ['commonArea', 'nonoverlappingArea', 'noncommonArea'], ascending = [False, True, True])
	return df



if __name__ == "__main__":
	pass
