import math
import multiprocessing
from typing import Dict, Generator, List, Optional, Tuple

import pandas
from loguru import logger
from tqdm import tqdm

try:
	from muller.clustering.metrics import distance_methods
	from muller import widgets
except ModuleNotFoundError:
	from . import distance_methods
	from ... import widgets

FilterType = Tuple[Optional[pandas.Series], Optional[pandas.Series]]


class DistanceCalculator:
	""" Refactor to clean up code and add multithreading."""

	def __init__(self, detection_limit: float, fixed_limit: float, metric: str, threads: Optional[int] = None):
		self.detection_limit = detection_limit
		self.fixed_limit = fixed_limit
		self.metric = metric
		self.threads = threads
		# Basically used as a cache. Should save memory compared to loading each pair of trajectories directly into `pair_combinations`.
		self.trajectories: Optional[pandas.DataFrame] = None

		self.progress_bar_minimum_points = 10000  # The value to activate the scale bar at.

	def calculate_pairwise_distances_threaded(self, pair_combinations: Generator, total: Optional[int] = None) -> Dict[
		Tuple[str, str], float]:
		""" Threaded version of the pairwise distance calculator"""
		pair_array: Dict[Tuple[str, str], float] = dict()
		progress_bar = tqdm(total = total)
		pool = multiprocessing.Pool(processes = self.threads)
		for i in tqdm(
				[pool.apply_async(calculate_distance, args = (self, e, self.trajectories)) for e in pair_combinations]):
			key, value = i.get()  # retrieve the calculated value
			pair_array[key] = value
			pair_array[key[::-1]] = value
			progress_bar.update(1)

		return pair_array

	def calculate_pairwise_distances_serial(self, pair_combinations: Generator, total: Optional[int] = None):
		""" Nonthreaded version of the pairwise distance calculator"""
		pair_array: Dict[Tuple[str, str], float] = dict()

		# The progressbar is not really useful if there aren;t a lot of combinations, since calculating the pairwise
		# distances is pretty fast. So disable the progressbar when the expected time to calculate all distances
		# Is less than, say, 20-30s.
		use_progressbar = total >= 10_000
		if use_progressbar:
			progress_bar = tqdm(total = total)

		for element in pair_combinations:
			key, value = calculate_distance(self, element, self.trajectories)
			pair_array[key] = value
			pair_array[key[1], key[
				0]] = value  # It's faster to add the reverse key rather than trying trying to get  test forward and reverse keys
			if use_progressbar:
				progress_bar.update(1)
		return pair_array

	def calculate_pairwise_distances(self, labels: List[str]) -> Dict[Tuple[str, str], float]:
		""" Implements the actual loop over all pairs of trajectories.
		"""
		# May as well move the combination function here so we don't have to pass an additional parameter specifying the total number
		# of trajectories so tqdm workd properly.

		total_elements = len(labels)
		total_combinations = widgets.calculate_number_of_combinations(total_elements)

		logger.debug(
			f"Generating pairwise combinations with {total_elements} items resulting in {total_combinations} combinations...")

		if total_combinations > 100_000_000:
			message = f"The provided dataset has {total_elements} trajectories, which requires {total_combinations} distance calculations." \
					  "This will require a long time to process (you may need to adjust the number of available threads with the --threads option)" \
					  "and will consume a large amount of memory (i.e. more than 16GB)."
			logger.warning(message)
		pair_combinations = widgets.get_pair_combinations(labels)

		if self.threads and self.threads > 1:  # One process is slower than using the serial method.
			logger.debug(f"Using multithreading...")
			pair_array = self.calculate_pairwise_distances_threaded(pair_combinations, total_combinations)
		else:
			logger.debug(f"Using a single thread...")
			pair_array = self.calculate_pairwise_distances_serial(pair_combinations, total_combinations)

		# Assume that any pair with NAN values are the maximum possible distance from each other.
		try:
			maximum_distance = max(filter(lambda s: not math.isnan(s), pair_array.values()))
		except ValueError:
			# Usually cause by max() being used on an empty sequence.
			message = f"Could not calculate the pairwise distances due to invalid series (usually because all measurements are below the detectionlimit"
			raise ValueError(message)

		pair_array = {k: (v if not math.isnan(v) else maximum_distance) for k, v in pair_array.items()}

		return pair_array

	def run(self, trajectories: pandas.DataFrame) -> Dict[Tuple[str, str], float]:
		"""
			Calculates the distance between all pairwise combinations of mutational trajectories. The total number of pairs can be calculated by
			`n!/(n-2)!`, which can be simlified to `n*(n-1)`.
			10 trajectories -> 90
			100 trajectories -> 9900
			1000 trajectories -> 999000
			31000 trajectories -> 9.61E8

		Parameters
		----------
		trajectories
		"""
		logger.debug("Calculating the pairwise values...")
		logger.debug(f"\t detection limit: {self.detection_limit}")
		logger.debug(f"\t fixed limit: {self.fixed_limit}")
		logger.debug(f"\t metric: {self.metric}")
		logger.debug(f"\t threads: {self.threads}")

		self.trajectories = trajectories

		pairwise_distances = self.calculate_pairwise_distances(trajectories.index)

		return pairwise_distances


def get_pair_category(left: pandas.Series, right: pandas.Series, dlimit: float, flimit: float) -> str:
	"""
		Categorizes the pair of mutational trajectories based on the dynamics of the
		measured frequencies over time.
		Each pair of mutational trajectories can be categorized into one of the following groups:
		- "onlyFixed": Both trajectories were only detected during timepoints where they were fixed.
		- "PartiallyFixed":Both trajectories were fixed, but only one trajectory was fixed during all timepoints where it was detected.
		- "bothFixed": Both trajectories where fixed during at least one sampled timepoint, and both had intermediate values.
		- "oneFixed": Only one of the trajectories was fixed.
		- "notFixed: Neither trajectory where fixed at any of the sampled timepoints.
	Returns
	-------
	str
		- 'onlyFixed': Both trajectories were only detected after they fixed.
		- 'partiallyFixed': Only one trajectory had intermediate values.
		- 'bothFixed': Both trajectories had non-fixed timepoints.
		- 'oneFixed': Only one trajectory ever fixed.
		- 'notFixed': Neither trajectory was ever fixed.
	"""
	fixed_left = widgets.get_fixed(left, flimit)
	fixed_right = widgets.get_fixed(right, flimit)

	left_was_fixed = len(fixed_left) > 0
	right_was_fixed = len(fixed_right)
	# Check if there are any intermediate values between the detection limit and the fixed limit.
	intermediate_left = widgets.get_intermediate(left, dlimit, flimit)
	intermediate_right = widgets.get_intermediate(right, dlimit, flimit)

	#undetected_left = widgets.get_undetected(left, dlimit)
	#undetected_right = widgets.get_undetected(right, dlimit)

	if left_was_fixed and right_was_fixed:
		# If both of these fixed at some point, the regions where they're both fixed means they should be
		# included within the same genotype. If, however, they both are both fixed but at separate timepoints,
		# They probably shouldn't be grouped together. It is also possible for one to be fixed and the other to
		# be only fixed.
		left_was_only_fixed = intermediate_left.empty
		right_was_only_fixed = intermediate_right.empty

		if left_was_only_fixed and right_was_only_fixed:
			# Both series only have undetected or fixed values
			category = 'onlyFixed'
		elif not left_was_only_fixed and not right_was_only_fixed:
			# Both series had intermediate values
			category = 'bothFixed'
		else:
			# One series only had fixed timepoints, the other did not.
			category = 'partiallyFixed'

		# So, check if the 'fixed' regions overlap.
		# left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, process.detection_limit, process.fixed_limit,
		#	inner = False)
		overlapping_regions = widgets.get_overlap_regions(left, right, flimit)
		# Need to make sure there aren't any regions after the first timepoint where both are fixed that is not fixed
		first_overlap_timepoint = widgets.get_first_fixed_timepoint(overlapping_regions, flimit)
		# Now get the section of the series after the first fixed timepoint for both series.
		overlapping_regions_after_fixed: pandas.Series = overlapping_regions[first_overlap_timepoint:]
		# Are there any timepoints after this where they aren't both overlapping?
		nonfixed_regions_after_first = overlapping_regions.apply(lambda s: s < 0.97)
		# If there is no subsequent timepoint that is not fixed or the overlapping regions after fixed is empty
		# Then these should be grouped together, but still need to check if the nonfixed regions agree.
		#is_genotype = overlapping_regions_after_fixed.empty or nonfixed_regions_after_first.sum() == 0

	# Check if the number of timepoints that overlap is more than one (overlapping suring one timepoint may be
	# due to measurement error.
	elif (left_was_fixed and not right_was_fixed) or (not left_was_fixed and right_was_fixed):
		# Only one trajectory had fixed timepoints
		category = 'oneFixed'
	elif not left_was_fixed and not right_was_fixed:
		# Neither series had fixed values.
		category = 'notFixed'
	else:
		message = f"Could not determine the category of this mutational pair"
		logger.error(message)
		logger.error(left)
		logger.error(right)
		raise ValueError

	return category


def filter_timepoints(left_trajectory: pandas.Series, right_trajectory: pandas.Series, dlimit: float,
		flimit: float) -> FilterType:
	"""
		Filters the available timepoints based on the measured dynamics.
	"""
	"""
	Legacy code:
	if left_was_fixed == right_was_fixed:
		left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, process.detection_limit, process.fixed_limit,
			inner = False)
	else:
		left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, process.detection_limit, inner = False)

	
	"""
	pair_category = get_pair_category(
		left_trajectory, right_trajectory,
		dlimit = dlimit, flimit = flimit
	)
	if pair_category == 'onlyFixed':
		left_reduced = right_reduced = None

	elif pair_category == 'partiallyFixed':
		# There is no overlap between these series so we have to rely on the overlap between "fixed" regions.
		left_reduced, right_reduced = widgets.get_valid_points(
			left_trajectory, right_trajectory,
			dlimit = dlimit, flimit = flimit, inner = False
		)

	elif pair_category == 'oneFixed':
		left_reduced, right_reduced = widgets.get_valid_points(
			left_trajectory, right_trajectory,
			dlimit = dlimit, inner = False
		)
	elif pair_category == 'notFixed':
		left_reduced, right_reduced = widgets.get_valid_points(
			left_trajectory, right_trajectory,
			dlimit = dlimit, flimit = flimit,
			inner = False
		)
	elif pair_category == 'bothFixed':
		left_reduced, right_reduced = widgets.get_valid_points(
			left_trajectory, right_trajectory,
			dlimit = dlimit, flimit=flimit,
			inner = False
		)
	else:
		message = f"Got an invalid category for a pair of trajectories: '{pair_category}'"
		logger.error(message)
		logger.error(left_trajectory.tolist())
		logger.error(right_trajectory.tolist())
		raise ValueError(message)
	return left_reduced, right_reduced


# Keep this as a separate function. Class methods are finicky when used with multiprocessing.
def calculate_distance(process: DistanceCalculator, element: Tuple[str, str], trajectories: pandas.DataFrame) -> Tuple[
	Tuple[str, str], float]:
	""" Implements the actual calculation for a specific pair of trajectories.
		It should be atomitized so that it works with multithreading.
	"""
	left, right = element
	left_trajectory = trajectories.loc[left]
	right_trajectory = trajectories.loc[right]

	# We only care about the timepoints such that `detection_cutoff` < f < `fixed_cutoff`.
	# For now, lets require that both timepoints are detected and not yet fixed.
	# There is an issue related to comparing fixed genotypes against non-fixed genotypes.

	left_reduced, right_reduced = filter_timepoints(
		left_trajectory, right_trajectory, process.detection_limit, process.fixed_limit
	)

	if left_reduced is None or right_reduced is None:
		# Treat both trajectories as fixed immediately.
		distance_between_series = fixed_overlap(left_trajectory, right_trajectory, process.fixed_limit)
	else:
		distance_between_series = distance_methods.calculate_distance(left_reduced, right_reduced, process.metric)

	return element, distance_between_series


def fixed_overlap(left: pandas.Series, right: pandas.Series, fixed_cutoff: float) -> float:
	"""
		Calculates the overlap of two series that are both undetected or fixed at all timepoints.
		This will be the distance value for these two series.
	Parameters
	----------
	left:pandas.Series
	right:pandas.Series
	fixed_cutoff: float
	"""
	left_fixed: pandas.Series = left[left.ge(fixed_cutoff)].dropna()
	right_fixed: pandas.Series = right[right.ge(fixed_cutoff)].dropna()
	if left_fixed.empty or right_fixed.empty:
		# Both are undetected, since they both failed the invalid timepoint filter.
		# As in, was already tested for intermediate values, so the only timepoints left would be undetected or fixed,
		# and these were just tested for fixed values. However, this situation is usually filtered during the earlier
		# steps so don't expect this if branch to be run.
		result = False
	else:
		# Since these genotypes fixed immediately, they should only be grouped together if they
		# fixed at the same timepoint.
		fixed_at_same_time = left_fixed.index[0] == right_fixed.index[0]
		# Check if the trajectories overlap completely
		overlap = set(left_fixed.index) & set(right_fixed.index)
		# This could probably be replaced with set(a) == set(b), but for some reason I compare them to
		# The union instead. It will do the same job and seems to work so I'll leave it for now.
		complete_overlap = len(overlap) == len(left_fixed) and len(overlap) == len(right_fixed)
		result = fixed_at_same_time and complete_overlap
	result = 0 if result else math.nan

	return result

