from typing import Any, List, Tuple

import pandas

try:
	from muller.muller_genotypes import generate
except ModuleNotFoundError:
	from muller_genotypes import generate


def get_fuzzy_backgrounds(genotypes: pandas.DataFrame, cutoffs: List[float]) -> Tuple[pandas.DataFrame, Tuple[float, float]]:
	""" Extracts the backgrounds using a list of frequency breakpoints."""
	for cutoff in cutoffs:
		backgrounds = genotypes[genotypes.max(axis = 1) > cutoff]
		backgrounds = backgrounds[backgrounds.sum(axis = 1) > 2] # Check if background appears at multiple timepoints.
		if not backgrounds.empty:
			fuzzy_fixed_cutoff = cutoff
			fuzzy_detected_cutoff = 1 - cutoff
			break
	else:
		raise ValueError("The filters cannot be applied since no backgrounds can be detected.")
	return backgrounds, (fuzzy_detected_cutoff, fuzzy_fixed_cutoff)


# noinspection PyTypeChecker
def check_if_genotype_is_invalid(genotype: pandas.Series, background_detected: int, background_fixed: int, detection_cutoff: float,
		use_strict_filter: bool) -> bool:
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

	if was_detected_before_and_after_background and was_detected_before_and_after_fixed:
		# To confirm that it is an invalid genotype rather than a genotype that was wiped out by a background and then reapeared,
		# Check to see if it was undetected at the timpont the background fixed.

		if use_strict_filter or genotype.loc[background_detected] > detection_cutoff:
			# The genotype was detected at the first fixed timepoint.
			return True
	return False


# noinspection PyTypeChecker
def get_first_timpoint(series: pandas.Series, cutoff: float) -> int:
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
def get_invalid_genotype(genotypes: pandas.DataFrame, detection_cutoff: float, cutoffs: List[float], use_strict_filter: bool) -> str:
	"""	Invalid genotypes are those that don't make sense in the context of evolved populations. For example, when a genotype fixes it wipes out all
		unrelated diversity and essentially 'resets' the mutation pool. Genotypes which are detected prior to a fixed genotype should, in theory,
		fall to an undetected frequency. Any genotypes that do not follow this rule (are detected both before and after a genotype fixes) should
		be considered invalid and removed from the population. The genotypes then should be re-calculated with the offending trajectories removed.
	Parameters
	----------
	genotypes: pands.DataFrame
	detection_cutoff: float
	cutoffs: List[float]
	use_strict_filter: bool
	Returns
	-------

	"""
	# Find all the backgrounds for this population. Some may fall below the usual `fixed_cutoff` threshold, so use the same frequency breakpoints
	# used when sorting the genotypes.
	backgrounds, (fuzzy_detected_cutoff, fuzzy_fixed_cutoff) = get_fuzzy_backgrounds(genotypes, cutoffs)
	# Make sure the selected detection cutoff is a valid float.
	fuzzy_detected_cutoff = max(detection_cutoff, fuzzy_detected_cutoff)

	# Filter out backgrounds that are only detected at one timepoint.
	backgrounds = _get_backgrounds_present_at_multiple_timepoints(backgrounds, fuzzy_detected_cutoff)

	# We want to iterate over the non-background genotypes to check if they appear both before and after any genotypes that fix.
	not_backgrounds = genotypes[~genotypes.index.isin(backgrounds.index)]

	# Iterate over the detected backgrounds.
	for _, background in backgrounds.iterrows():
		# Find the timepoint where the background first fixes.
		first_detected_point: int = get_first_timpoint(background, fuzzy_detected_cutoff)
		first_fixed_point: int = get_first_timpoint(background, fuzzy_fixed_cutoff)

		# Iterate over the non-background genotypes.
		for genotype_label, genotype in not_backgrounds.iterrows():
			# Double check that it is not a background
			if genotype_label in backgrounds.index: continue
			# Check if it is invalid.
			is_invalid = check_if_genotype_is_invalid(genotype, first_detected_point, first_fixed_point, fuzzy_detected_cutoff, use_strict_filter)
			if is_invalid:
				return genotype_label


DF = pandas.DataFrame


def filter_genotypes(trajectory_table: pandas.DataFrame, goptions: generate.GenotypeOptions, frequency_cutoffs: List[float],
		use_strict_filter: bool) -> Tuple[DF, DF, Any, Any]:
	"""
		Iteratively calculates the population genotypes, checks and removes invalid genotypes, and recomputes the genotypes until no changes occur.
	Parameters
	----------
	trajectory_table: pandas.DataFrame
	goptions: GenotypeOptions
	frequency_cutoffs: List[float]
	use_strict_filter: bool

	Returns
	-------

	"""
	# Remove 1.0 fro mthe list of frequency breakpoints to account for measurement errors.
	frequency_cutoffs = [i for i in frequency_cutoffs if i <= goptions.fixed_breakpoint]
	trajectory_table = trajectory_table.copy(deep = True)  # To avoid any unintended changes to the original table.
	genotype_table, genotype_members, linkage_table = generate.generate_genotypes(trajectory_table, options = goptions)

	# cache: List[Tuple[DF, DF]] = [(trajectory_table.copy(), genotype_table.copy())]
	_iterations = 20  # arbitrary, used to ensure the program does not encounter an infinite loop.
	for _ in range(_iterations):
		# Search for genotypes that do not make sense in the context of an evolved population.
		current_invalid_genotype = get_invalid_genotype(genotype_table, goptions.detection_breakpoint, frequency_cutoffs, use_strict_filter)
		if current_invalid_genotype is None:
			break
		else:
			# Get a list of the trajectories that form this genotype.
			invalid_members = genotype_members.loc[current_invalid_genotype].split('|')
			# Remove these trajectories from the trajectories table.
			trajectory_table = trajectory_table[~trajectory_table.index.isin(invalid_members)]
			# Re-calculate the genotypes based on the remaining trajectories.
			genotype_table, genotype_members, linkage_table = generate.generate_genotypes(trajectory_table, options = goptions)

	# cache.append((trajectory_table.copy(), genotype_table.copy()))
	# Update the trajectories that comprise each genotype.
	else:
		print(f"Could not filter the genotypes after {_iterations} iterations.")
	return trajectory_table, genotype_table, genotype_members, linkage_table


if __name__ == "__main__":
	pass
