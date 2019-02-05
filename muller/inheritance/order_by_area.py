import pandas
from scipy.integrate import simps
import numpy
from functools import partial
from typing import Set
def area_of_series(series: pandas.Series) -> float:
	""" Calculates the area of a discrete series."""
	#area = numpy.trapz(series.values)#, series.index)
	#area = simps(series.values, series.index)
	#area = series.sum()
	total = series.sum()
	return total


def calculate_noncommon_area(left: pandas.Series, right: pandas.Series, dlimit:float) -> float:
	""" Calculates the area not contained by bboth `left` and `right`"""
	difference = abs(left - right)
	noncommon_area = area_of_series(difference)

	return noncommon_area

def calculate_non_overlapping_area(left:pandas.Series, right:pandas.Series, dlimit)->float:
	# Area of the timepoionts that do no overlap
	overlap_index = get_overlap(left, right, dlimit)
	nonoverlap_left = left.loc[~left.index.isin(overlap_index)]
	nonoverlap_right = right.loc[~right.index.isin(overlap_index)]
	area_left = area_of_series(nonoverlap_left)
	area_right = area_of_series(nonoverlap_right)
	return area_left + area_right
def calculate_common_area(left: pandas.Series, right: pandas.Series, dlimit: float) -> float:
	#left = left.reset_index(drop = True)
	#right = right.reset_index(drop = True)
	overlap_index = get_overlap(left, right, dlimit)
	if not overlap_index:
		overlap = 0
	else:
		overlap_left_series = left.loc[overlap_index]
		overlap_right_series = right.loc[overlap_index]
		overlap_left = area_of_series(overlap_left_series)
		overlap_right = area_of_series(overlap_right_series)

		overlap = area_of_series(left.where(left < right, right))
	return overlap


def get_overlap(left: pandas.Series, right: pandas.Series, dlimit: float) -> Set[int]:
	detected_left = left[left > dlimit]
	detected_right = right[right > dlimit]

	overlap_index = set(detected_left.index) & set(detected_right.index)
	return overlap_index

def compare_to_table(left:pandas.Series, table:pandas.DataFrame, dlimit:float)->pandas.DataFrame:
	area = area_of_series(left)
	axis = 1
	common_area = table.apply(calculate_common_area, axis = axis, args = (left, dlimit))/ area
	noncommon_area = table.apply(calculate_noncommon_area, axis = axis, args = (left, dlimit))/ area
	nonoverlapping_area = table.apply(calculate_non_overlapping_area, axis = axis, args = (left, dlimit))/ area

	df = pandas.concat([common_area, noncommon_area, nonoverlapping_area], axis = 1)
	df.columns = ['commonArea', 'noncommonArea', 'nonoverlappingArea']
	df = df.sort_values(by = ['commonArea', 'nonoverlappingArea', 'noncommonArea'], ascending = [False, True, True])
	return df



if __name__ == "__main__":
	from import_data import import_table_from_string
	string = """
		Genotype	0	17	25	44	66	75	90
		genotype-1	0	0	0.261	1	1	1	1
		genotype-6	0	0	0	0.273	0.781	1	1
		genotype-3	0	0	0	0	0	1	1
		genotype-4	0	0	0	0.525	0.454	0.911	0.91
		genotype-5	0	0	0	0.147	0.45	0.924	0.887
		genotype-11	0	0	0	0	0.278	0.822	0.803
		genotype-2	0	0.38	0.432	0	0	0	0
		genotype-8	0	0	0	0.403	0.489	0.057	0.08
		genotype-14	0	0	0	0	0	0.2675	0.326
		genotype-10	0	0	0	0.138	0.295	0	0.081
		genotype-12	0	0	0	0	0.2335	0.133	0.0375
		genotype-7	0	0	0	0.188	0.171	0.232	0.244
		genotype-9	0	0	0.117	0	0	0	0.103
		genotype-13	0	0	0.033	0.106	0.1065	0	0
		genotype-15	0	0	0	0.1145	0	0.1205	0.0615
	"""
	table = import_table_from_string(string, index = 'Genotype')
	table.columns = [int(i) for i in table.columns]
	left = table.loc['genotype-2']

	result = compare_to_table(left, table, 0.03)
	print(result.to_string())

