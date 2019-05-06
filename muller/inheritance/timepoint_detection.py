import pandas


# noinspection PyTypeChecker,PyUnresolvedReferences
def _get_timepoint_above_threshold(transposed_timepoints: pandas.DataFrame, cutoff: float, name: str = None) -> pandas.Series:
	"""
		Calculates when a genotype was first fixed based on `cutoff`. The resulting series will be named `name`, if given.
	Parameters
	----------
	transposed_timepoints: pandas.DataFrame
		A dataframe indexed by timepoint and with each series a separate column.
	cutoff: float
		The threshold a frequency must exceed to be considered a valid timepoint.

	Returns
	-------
	pandas.Series
		A series mapping series name to the first timepoint that series first exceeded the threshold value.
		Series that never exceed the threshold are mapped to a value of 0.
	"""
	threshold_series: pandas.Series = (transposed_timepoints > cutoff).idxmax(0).sort_values()
	if name:
		threshold_series.name = name
	return threshold_series


def get_first_detected_timepoint(transposed_genotypes: pandas.DataFrame, cutoff: float) -> pandas.Series:
	"""

	Parameters
	----------
	transposed_genotypes:pandsa.DataFrame
		A table where rows correspond to timepoints and columns correspond to genotypes.
	cutoff

	Returns
	-------

	"""
	initial_genotype_values = transposed_genotypes.iloc[0].transpose()
	first_detected = _get_timepoint_above_threshold(transposed_genotypes, cutoff, 'firstDetected')
	first_detected_reduced = first_detected[(first_detected != 0) | (initial_genotype_values > cutoff)]

	return first_detected_reduced


def get_first_significant_timepoint(transposed_genotypes: pandas.DataFrame, cutoff: float) -> pandas.Series:
	"""
		Retrieves the first timepoint in each genotye that exceeds the significant cutoff.
	Parameters
	----------
	transposed_genotypes
	cutoff

	Returns
	-------
	"""

	initial_genotype_values = transposed_genotypes.iloc[0].transpose()
	first_above_threshold = _get_timepoint_above_threshold(transposed_genotypes, cutoff, 'firstSignificant')
	first_above_threshold_reduced_dict = dict()
	for key, value in first_above_threshold.items():
		initial_value = initial_genotype_values[key]
		if value == transposed_genotypes.index[0]:
			value = value if initial_value > cutoff else transposed_genotypes.index[-1]
		first_above_threshold_reduced_dict[key] = value

	first_above_threshold_reduced = pandas.Series(first_above_threshold_reduced_dict, name = 'firstThreshold')
	return first_above_threshold_reduced


def get_first_fixed_timepoint(transposed_genotypes: pandas.DataFrame, cutoff: float) -> pandas.DataFrame:
	"""
		Retrives the first timepoint every genotype exceeded the cutoff value.
	Parameters
	----------
	transposed_genotypes
	cutoff

	Returns
	-------

	"""
	first_fixed = _get_timepoint_above_threshold(transposed_genotypes, cutoff, 'firstFixed')
	if cutoff != 0:
		# If we are checking for the lowest frequency breakpoint, ignore this check.
		first_fixed_reduced: pandas.DataFrame = first_fixed.iloc[first_fixed.to_numpy().nonzero()]
	else:
		first_fixed_reduced = first_fixed.to_frame()

	return first_fixed_reduced
