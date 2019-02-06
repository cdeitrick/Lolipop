from typing import Any, Optional, Tuple

import pandas

try:
	from inheritance import order_by_area
except ModuleNotFoundError:
	from . import order_by_area

def get_detected_points(left: pandas.Series, right: pandas.Series, dlimit: float, exclude_fixed:bool = False) -> pandas.DataFrame:
	""" Removes points where at least one series was not detected."""
	df = pandas.concat([left, right], axis = 1)
	df.columns = ['left', 'right'] # Prevent an IndexError when `left` and `right` refer to the same series.
	at_least_one_detected = df.sum(axis = 1) > dlimit
	at_least_one_detected = at_least_one_detected[at_least_one_detected]
	if len(at_least_one_detected) == 1:
		result = df.loc[at_least_one_detected.index[0]]
	else:
		result = df.loc[at_least_one_detected.index]

	if exclude_fixed:
		not_fixed = (result['left'] < 0.97) & (result['right'] < 0.97)
		result = result.loc[not_fixed[not_fixed].index]

	return result


def check_additive_background(left: pandas.Series, right: pandas.Series, double_cutoff: float, single_cutoff: float) -> bool:
	"""
		Want to check if the sum of two genotypes at one or more timepoints is consistently greater than 1.0
	Parameters
	----------
	left, right: pandas.Series
	double_cutoff: float
		Returns True if two or more timepoints sum to a value that exceeds this frequency.
	single_cutoff
		Returns True if one or more timepoints sum to a value that exceeds this frequency.
	Returns
	-------
	bool
	"""
	#trajectorysum = right + left
	#trajectorysum = trajectorysum[trajectorysum > 0.03]
	trajectorysum = get_detected_points(left, right, 0.03).sum(axis = 1)
	double_check = (trajectorysum > double_cutoff).sum() > 0  # Implicit conversion from bool to int.
	single_check = (trajectorysum > single_cutoff).sum() > 0
	return double_check or single_check


def check_subtractive_background(left: pandas.Series, right: pandas.Series, double_cutoff: float,
		single_cutoff: float) -> bool:
	"""
		Check if the one genotype is consistently and significantly larger than the other.
	Parameters
	----------
	left: pandas.Series
	right: pandas.Series
	double_cutoff: float
	single_cutoff: float

	Returns
	-------

	"""
	# Check if the current genotype is over 15% larger than the other.
	diff_trajectory = -abs(right - left)
	double_diff_trajectory = (diff_trajectory < double_cutoff).sum() > 1  # At least two timepoints exceed left.
	single_diff_trajectory = (diff_trajectory < single_cutoff).sum() > 0  # implicit conversion from bool to int
	return double_diff_trajectory or single_diff_trajectory


# noinspection PyTypeChecker
def check_derivative_background_legacy(left: pandas.Series, right: pandas.Series, detection_cutoff: float) -> float:
	# Look at the first point when both trajectories are non-zero, and observe delta over time.
	# Are derivatives correlated or anti-correlated? If sufficiently correlated, they are on the same background.
	right_trajectory_filtered: pandas.Series = right > detection_cutoff  # All positions > 0.02 will evaluate to True
	left_trajectory_filtered: pandas.Series = left > detection_cutoff
	# Find the points where neither left nor right is zero.
	startpoint = max([right_trajectory_filtered.idxmax(), left_trajectory_filtered.idxmax()])
	endpoint = min([right_trajectory_filtered[::-1].idxmax(), left_trajectory_filtered[::-1].idxmax()])

	# Convert the index labels to their corresponding integer position in the index.
	# Pandas can't use labels to get a slice of a series.
	start_index = left.index.get_loc(startpoint)
	end_index = left.index.get_loc(endpoint)
	# The matlab script uses the difference of the following timepoint versus the curent one,
	# while pandas default behaviour is current timepoint - previous timepoint.
	# Need to indicate difference of following timepoint, and reverse polarity(negative sign)
	right_derivative = right.diff(-1)
	left_derivative = left.diff(-1)

	right_window = right_derivative[start_index:end_index] * -1
	left_window = left_derivative[start_index:end_index] * -1
	delta = left_window.dot(right_window)

	return delta


def check_derivative_background(left: pandas.Series, right: pandas.Series, detection_cutoff: float) -> float:
	# Pandas implementation of the derivative check, since it basically just checks for covariance.
	df = get_detected_points(left, right, detection_cutoff, exclude_fixed = True)
	#return left.cov(right)
	return df.iloc[:, 0].cov(df.iloc[:, 1])

def apply_genotype_checks(type_trajectory: pandas.Series, test_trajectory: pandas.Series, options: Any) -> Tuple[bool, bool, Optional[float]]:
	""" Applies the three checks to `type_trajectory` and `test_trajectory`."""
	additive_check = check_additive_background(
		left = type_trajectory,
		right = test_trajectory,
		double_cutoff = options.additive_background_double_cutoff,
		single_cutoff = options.additive_background_single_cutoff
	)

	subtractive_check = check_subtractive_background(
		left = type_trajectory,
		right = test_trajectory,
		double_cutoff = options.subtractive_background_double_cutoff,
		single_cutoff = options.subtractive_background_single_cutoff
	)

	if subtractive_check:
		delta = None
	else:
		delta = check_derivative_background(
			left = type_trajectory,
			right = test_trajectory,
			detection_cutoff = options.derivative_detection_cutoff
		)

	return additive_check, subtractive_check, delta


def apply_genotype_checks_to_table(unnested_trajectory: pandas.Series, table: pandas.DataFrame, options: Any) -> pandas.DataFrame:
	axis = 1
	additive_check = table.apply(
		check_additive_background,
		args = (unnested_trajectory, options.additive_background_double_cutoff, options.additive_background_single_cutoff),
		axis = axis
	)

	subtractive_check = table.apply(
		check_subtractive_background,
		args = (unnested_trajectory, options.subtractive_background_double_cutoff, options.subtractive_background_single_cutoff),
		axis = axis
	)
	covariance_check = table.apply(
		check_derivative_background,
		args = (unnested_trajectory, options.derivative_detection_cutoff),
		axis = axis
	)

	area_check = table.apply(
		order_by_area.calculate_common_area,
		args = (unnested_trajectory, 0.03),
		axis = axis
	)
	area_check = area_check / order_by_area.area_of_series(unnested_trajectory)

	columns = ['additiveCheck', 'subtractiveCheck', 'derivativeCheck', 'areaCheck']

	df = pandas.concat([additive_check, subtractive_check, covariance_check, area_check], axis = 1)
	df.columns = columns
	return df


if __name__ == "__main__":
	from import_data import import_table_from_string
	import argparse

	args = argparse.Namespace()
	args.additive_background_double_cutoff = 1.03
	args.additive_background_single_cutoff = 1.15
	args.subtractive_background_double_cutoff = -0.03
	args.subtractive_background_single_cutoff = -0.15
	args.derivative_detection_cutoff = 0.01

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
	left = table.loc['genotype-4']
	right = table.loc['genotype-3']

	df = apply_genotype_checks_to_table(left, table, args)
	print(df.to_string())
	print(left.cov(right))
	print(check_derivative_background_legacy(left, right, 0.03))
