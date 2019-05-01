
from typing import List, Tuple

import pandas

from loguru import logger
try:
	from muller.options import GenotypeOptions
except ModuleNotFoundError:
	from options import GenotypeOptions


def get_fuzzy_backgrounds(genotypes: pandas.DataFrame, cutoffs: List[float]) -> Tuple[pandas.DataFrame, float]:
	""" Extracts the backgrounds using a list of frequency breakpoints."""
	for cutoff in cutoffs:
		backgrounds = genotypes[genotypes.max(axis = 1) > cutoff]
		backgrounds = backgrounds[backgrounds.sum(axis = 1) > 2]  # Check if background appears at multiple timepoints.
		if not backgrounds.empty:
			fuzzy_fixed_cutoff = cutoff
			break
	else:
		raise ValueError("The filters cannot be applied since no backgrounds can be detected.")
	return backgrounds, fuzzy_fixed_cutoff


# noinspection PyTypeChecker
def check_if_genotype_is_invalid(genotype: pandas.Series, background_detected: int, background_fixed: int, detection_cutoff: float,
		background_value: float, use_strict_filter: bool) -> bool:
	"""
		Checks if a genotype does not adhere to certain assumptions related to the background.
		1. The genotype should not be present both before and after a background fixes, and should be nonzero when the background fixes.
	Parameters
	----------
	genotype: pandas.Series
		The genotype being tested.
	background_detected: int
		The first timepoint the background was detected.
	background_fixed: int
		The first timepoint the background qualifies as 'fixed'
	detection_cutoff: float
		The cutoff to determine whether a trajectory counts as 'detected'. It is based on the frequency cutoff used to identify the backgrounds.
	background_value: float
		The value of the background when it is considered fixed. This is needed since we are using fuzzy cutoffs to find the most abundant genotype.
	use_strict_filter: bool
		Whether to eliminate every genotype detected before and after a fixed point, or to allow these genotypes if they have frequency 0 at the fixed timepoint.

	Returns
	-------
	bool
	"""
	# get the nonzero timepoints for the genotype.
	detected_points = genotype[genotype > detection_cutoff]

	if len(detected_points) < 2:
		# The genotype is either undetected at all timepoints or is detected once.
		return False

	first_detected = detected_points.first_valid_index()
	last_detected = detected_points.last_valid_index()

	# Check if the genotype was detected prior to the background. Skip those that were not.
	was_detected_before_and_after_background = first_detected < background_detected < last_detected
	# Check if the genotype was detected before and after the timepoint the current backgound fixed.
	was_detected_before_and_after_fixed = first_detected < background_fixed < last_detected
	# print(genotype.name, (first_detected, background_detected), (last_detected, background_detected), was_detected_before_and_after_fixed, was_detected_before_and_after_background)
	if was_detected_before_and_after_background and was_detected_before_and_after_fixed:
		# To confirm that it is an invalid genotype rather than a genotype that was wiped out by a background and then reapeared,
		# Check to see if it was undetected at the timpont the background fixed.
		value_at_fixed_point = genotype.loc[background_fixed]
		fixed_point_value = value_at_fixed_point + background_value
		if use_strict_filter or (value_at_fixed_point > detection_cutoff and fixed_point_value > (1 + detection_cutoff)):
			# The genotype was detected at the first fixed timepoint.
			return True
	return False


# noinspection PyTypeChecker
def get_first_timepoint_above_cutoff(series: pandas.Series, cutoff: float) -> int:
	""" Extracts the first timepoint that exceeds the cutoff."""
	return series[series > cutoff].idxmin()


def _get_backgrounds_present_at_multiple_timepoints(backgrounds: pandas.DataFrame, detection_cutoff: float) -> pandas.DataFrame:
	""" Filters out backgrounds that only appeared at one timepoint.
		We want to exclude 'backgrounds' which are only present at one timepoint. This is more likely to be a measurement error.
		Since the fixed threshold is almost 1 (ex. 0.97), a simple test to check for detection at multiple timepoints is to see if the sum
		of the background frequencies is greater than 1."""
	fuzzy_backgrounds = backgrounds.sum(axis = 1) > (1 + detection_cutoff)
	backgrounds = backgrounds[fuzzy_backgrounds]
	return backgrounds


# noinspection PyTypeChecker
def find_first_invalid_genotype(genotypes: pandas.DataFrame, backgrounds: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float,
		use_strict_filter: bool) -> str:
	"""	Invalid genotypes are those that don't make sense in the context of evolved populations. For example, when a genotype fixes it wipes out all
		unrelated diversity and essentially 'resets' the mutation pool. Genotypes which are detected prior to a fixed genotype should, in theory,
		fall to an undetected frequency. Any genotypes that do not follow this rule (are detected both before and after a genotype fixes) should
		be considered invalid and removed from the population. The genotypes then should be re-calculated with the offending trajectories removed.
	Parameters
	----------
	genotypes: pands.DataFrame
	backgrounds: pandas.DataFrame
	detection_cutoff: float
	fixed_cutoff: float
	use_strict_filter: bool
	Returns
	-------

	"""

	# Filter out backgrounds that are only detected at one timepoint.
	# backgrounds = _get_backgrounds_present_at_multiple_timepoints(backgrounds, detection_cutoff)
	# We want to iterate over the non-background genotypes to check if they appear both before and after any genotypes that fix.
	not_backgrounds = genotypes[~genotypes.index.isin(backgrounds.index)]
	# Iterate over the detected backgrounds.
	for _, background in backgrounds.iterrows():
		# Find the timepoint where the background first fixes.
		first_detected_point: int = get_first_timepoint_above_cutoff(background, detection_cutoff)
		first_fixed_point: int = get_first_timepoint_above_cutoff(background, fixed_cutoff)
		background_value_at_first_fixed_point = background.loc[first_fixed_point]
		# Iterate over the non-background genotypes.
		for genotype_label, genotype in not_backgrounds.iterrows():
			# Double check that it is not a background
			if genotype_label in backgrounds.index: continue
			# Check if it is invalid.
			is_invalid = check_if_genotype_is_invalid(genotype, first_detected_point, first_fixed_point, detection_cutoff,
				background_value_at_first_fixed_point, use_strict_filter)
			if is_invalid:
				return genotype_label


def _remove_single_point_background(genotypes: pandas.DataFrame, dlimit: float, flimit: float) -> List[str]:
	"""
		Identifies backgrounds which are only detected at a only single timepoint. Returns a list of labels that fail the filter.
	Parameters
	----------
	genotypes: pandas.DataFrame
	dlimit: float
		The detection limit
	flimit: float
		The fixed detection limit.

	Returns
	-------
	pandas.DataFrame
		A list of backgrounds that fail the filter.
	"""
	backgrounds = genotypes[genotypes.max(axis = 1) > flimit]
	invalid_backgrounds = backgrounds[backgrounds.sum(axis = 1) < (flimit + (2 * dlimit))]

	return invalid_backgrounds.index


def filter_trajectories(trajectory_table: pandas.DataFrame, dlimit: float, flimit: float) -> pandas.DataFrame:
	"""
		Filters out individual trajectories that fail certain filters.
	Parameters
	----------
	trajectory_table:pandas.DataFrame
	dlimit: float
		The detection limit
	flimit: float
		The value above which a genotype is considered 'fixed'

	Returns
	-------

	"""
	# Remove trajectories that only exist at one timepoint and exceed the fixed cutoff limit.
	failed_single_point_test = _remove_single_point_background(trajectory_table, dlimit, flimit)
	logger.info("These trajectories did not pass the trajectory filters: " + str(list(failed_single_point_test)))
	return trajectory_table[~trajectory_table.index.isin(failed_single_point_test)]


def filter_genotypes(original_genotypes: pandas.DataFrame, genotype_members: pandas.Series,detection_breakpoint:float, breakpoints: List[float],
		use_strict_filter: bool = False) -> List[str]:
	"""
		Finds all trajectories which should be filtered out of the dataset based on certain criteria.
		 These criteria concern both the trajectories and their parent genotypes.
	Parameters
	----------
	original_genotypes: pandas.DataFrame
		pre-computed genotypes table.
	genotype_members: pandas.Series
		Maps genotypes to a '|' delimited string of member trajectories.
	detection_breakpoint: float
	breakpoints: List[float]
	use_strict_filter: bool; default False

	Returns
	-------
	List[str]
		A list of all trajectories which should be filtered out of the dataset.
	"""
	# Find all the backgrounds for this population. Some may fall below the usual `fixed_cutoff` threshold, so use the same frequency breakpoints
	# used when sorting the genotypes.
	try:
		current_backgrounds, fuzzy_fixed_limit = get_fuzzy_backgrounds(original_genotypes, breakpoints)
	except ValueError:
		return []

	logger.debug(f"Backgrounds for filtering:" + str(list(current_backgrounds.index)))

	# Search for genotypes that do not make sense in the context of an evolved population.
	current_invalid_genotype = find_first_invalid_genotype(
		original_genotypes,
		current_backgrounds,
		detection_breakpoint,
		fuzzy_fixed_limit,
		use_strict_filter
	)


	if current_invalid_genotype is None:
		return []
	else:
		logger.info(f"Found Invalid genotype: {current_invalid_genotype}. Removing...")
		# Get a list of the trajectories that form this genotype.
		invalid_members = genotype_members.loc[current_invalid_genotype].split('|')
		logger.info("Invalid members: " + str(invalid_members))
		return invalid_members


if __name__ == "__main__":
	pass
