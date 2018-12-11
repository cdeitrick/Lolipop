from pathlib import Path
from typing import Any, List, Tuple

import pandas

try:
	from muller.muller_genotypes import calculate_genotypes
except ModuleNotFoundError:
	from muller_genotypes import calculate_genotypes

pandas.set_option('display.width', 300)


# noinspection PyTypeChecker
def get_invalid_genotype(genotypes: pandas.DataFrame, detection_cutoff: float, fixed_cutoff: float, cutoffs: List[float]) -> pandas.DataFrame:
	"""	Invalid genotypes are those that don't make sense in the context of evolved populations. For example, when a genotype fixes it wipes out all
		unrelated diversity and essentially 'resets' the mutation pool. Genotypes which are detected prior to a fixed genotype should, in theory,
		fall to an undetected frequency. Any genotypes that do not follow this rule (are detected both before and after a genotype fixes) should
		be considered invalid and removed from the population. The genotypes then should be re-calculated with the offending trajectories removed.
	Parameters
	----------
	genotypes: pands.DataFrame
	detection_cutoff: float
	fixed_cutoff: float
	cutoffs: List[float]

	Returns
	-------

	"""
	# fuzzy_fixed_cutoff = 0.5
	# cutoffs = [fixed_cutoff] + [1, .9, .8, .7, .6, .5, .4, .3, .2, .1, 0.0]
	# Find all the backgrounds for this population. Some may fall below the usual `fixed_cutoff` threshold, so use the same frequency breakpoints
	# used when sorting the genotypes.
	for cutoff in cutoffs:
		backgrounds = genotypes[genotypes.max(axis = 1) > cutoff]
		if not backgrounds.empty:
			fuzzy_fixed_cutoff = cutoff
			break
	else:
		raise ValueError("The filters cannot be applied since no backgrounds can be detected.")

	# We want to exclude 'backgrounds' which are only present at one timepoint. This is more likely to be a measurement error.
	# Since the fixed threshold is almost 1 (ex. 0.97), a simple test to check for detection at multiple timepoints is to see if the sum
	# of the background frequencies is greater than 1.
	fuzzy_backgrounds = backgrounds.sum(axis = 1) > (1 + detection_cutoff)
	backgrounds = backgrounds[fuzzy_backgrounds]

	# Extract the timepoints each background first fixed at.
	fixed_timepoints: pandas.DataFrame = backgrounds[backgrounds > fuzzy_fixed_cutoff].dropna(how = 'all').transpose()
	fixed_timepoints = fixed_timepoints.dropna(how = 'all').transpose()
	fixed_timepoints = [s.first_valid_index() for _, s in fixed_timepoints.iterrows()]
	# We want to iterate over the non-background genotypes to check if they appear both before and after any genotypes that fix.
	not_backgrounds = genotypes[~genotypes.index.isin(backgrounds.index)]

	# Iterate over the detected backgrounds.
	for _, background in backgrounds.iterrows():
		# Find the timepoint where the background first fixes.
		fuzzy_fixed_timepoints = background[background > fuzzy_fixed_cutoff]
		fixed_point: int = fuzzy_fixed_timepoints.first_valid_index()

		# Iterate over the non-background genotypes.
		for genotype_label, genotype in not_backgrounds.iterrows():
			# Double check that it is not a background
			if genotype_label in backgrounds.index: continue

			# get the nonzero timepoints for the genotype.
			detected_points = genotype[genotype > detection_cutoff]

			if len(detected_points) < 2:
				# The genotype is either undetected at all timepoints or is detected once.
				continue

			first_detected = detected_points.first_valid_index()
			last_detected = detected_points.last_valid_index()

			# Check if the genotype was detected before and after the timepoint the current backgound fixed.
			if first_detected < fixed_point < last_detected:
				# To confirm that it is an invalid genotype rather than a genotype that was wiped out by a background and then reapeared,
				# Check to see if it was undetected at the timpont the background fixed.
				first_fixed_points: pandas.Series = genotype[fixed_timepoints]
				zero_points = first_fixed_points[first_fixed_points <= detection_cutoff]
				# The genotype was seen both before and after the background fixed.
				# Check if it was undetected when a genotype fixed.
				if zero_points.empty:
					return genotype_label
				else:
					# Remove even if undetected at fixed points.
					return genotype_label


DF = pandas.DataFrame


def workflow(trajectory_table: pandas.DataFrame, goptions: calculate_genotypes.GenotypeOptions, frequency_cutoffs: List[float]) -> Tuple[DF, DF, Any]:
	"""
		Iteratively calculates the population genotypes, checks and removes invalid genotypes, and recomputes the genotypes until no changes occur.
	Parameters
	----------
	trajectory_table: pandas.DataFrame
	goptions: GenotypeOptions
	frequency_cutoffs: List[float]

	Returns
	-------

	"""
	trajectory_table = trajectory_table.copy(deep = True) # To avoid any unintended changes to the original table.
	genotype_table = calculate_genotypes.workflow(trajectory_table, options = goptions)

	cache: List[Tuple[DF, DF]] = [(trajectory_table.copy(), genotype_table.copy())]
	original_genotype_members = genotype_table.pop('members')
	_iterations = 10 # arbitrary, used to ensure the program does not encounter an infinite loop.
	for _ in range(_iterations):
		# Search for genotypes that do not make sense in the context of an evolved population.
		current_invalid_genotype = get_invalid_genotype(genotype_table, goptions.detection_breakpoint, goptions.fixed_breakpoint, frequency_cutoffs)
		if current_invalid_genotype is None:
			break
		else:
			# Get a list of the trajectories that form this genotype.
			invalid_members = original_genotype_members.loc[current_invalid_genotype].split('|')
			# Remove these trajectories from the trajectories table.
			trajectory_table = trajectory_table[~trajectory_table.index.isin(invalid_members)]
			# Re-calculate the genotypes based on the remaining trajectories.
			genotype_table = calculate_genotypes.workflow(trajectory_table, options = goptions)

			cache.append((trajectory_table.copy(), genotype_table.copy()))
			# Update the trajectories that comprise each genotype.
			original_genotype_members = genotype_table.pop('members')
	else:
		print(f"Could not filter the genotypes after {_iterations} iterations.")
	genotype_table['members'] = original_genotype_members
	return trajectory_table, genotype_table, cache


if __name__ == "__main__":
	pass
