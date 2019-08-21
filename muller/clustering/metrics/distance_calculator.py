import itertools
import math
import multiprocessing
from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger

try:
	# tqdm is an optional dependancy
	from tqdm import tqdm
except ModuleNotFoundError:
	tqdm = None

try:
	from muller.clustering.metrics import distance_methods
	from muller import widgets
except ModuleNotFoundError:
	from . import distance_methods
	from ... import widgets


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

	def run(self, trajectories: pandas.DataFrame) -> Dict[Tuple[str, str], float]:
		logger.debug("Calculating the pairwise values...")
		logger.debug(f"\t detection limit: {self.detection_limit}")
		logger.debug(f"\t fixed limit: {self.fixed_limit}")
		logger.debug(f"\t metric: {self.metric}")
		logger.debug(f"\t threads: {self.threads}")

		self.trajectories = trajectories
		# Use a list so that we can get the size for tqdm. Also prevent weird errors if we use the variable after consuming it.
		# noinspection PyTypeChecker
		pair_combinations: List[Tuple[str, str]] = list(itertools.combinations(trajectories.index, 2))

		pairwise_distances = self.calculate_pairwise_distances(pair_combinations)

		return pairwise_distances

	def calculate_pairwise_distances(self, pair_combinations: List[Tuple[str, str]]) -> Dict[Tuple[str, str], float]:
		""" Implements the actual loop over all pairs of trajectoies. """

		# Initialize a progressbar if the dataset has a lot of combinations.
		# if len(pair_combinations) >= self.progress_bar_minimum_points and tqdm:
		#	self.progress_bar = tqdm(total = len(pair_combinations))
		# else:
		#	self.progress_bar = None

		pair_array: Dict[Tuple[str, str], float] = dict()

		if self.threads and self.threads > 1: #One process is slower than using the serial method.
			pool = multiprocessing.Pool(processes = self.threads)
			for i in tqdm([pool.apply_async(calculate_distance, args = (self, e, self.trajectories)) for e in pair_combinations]):
				key, value = i.get()
				pair_array[key] = value
		else:
			for element in pair_combinations:
				key, value = calculate_distance(self, element, self.trajectories)
				pair_array[key] = value
				pair_array[
					key[1], key[0]] = value  # It's faster to add the reverse key rather than trying trying to get  test forward and reverse keys

		# Assume that any pair with NAN values are the maximum possible distance from each other.
		maximum_distance = max(filter(lambda s: not math.isnan(s), pair_array.values()))
		pair_array = {k: (v if not math.isnan(v) else maximum_distance) for k, v in pair_array.items()}

		return pair_array


# Keep this as a separate function. Class methods are finicky when used with multiprocessing.
def calculate_distance(process, element: Tuple[str, str], trajectories) -> Tuple[Tuple[str, str], float]:
	""" Implements the actual calculation for a specific pair of trajectories.
		It should be atomitized so that it works with multithreading.
	"""
	left, right = element
	left_trajectory = trajectories.loc[left]
	right_trajectory = trajectories.loc[right]

	# We only care about the timepoints such that `detection_cutoff` < f < `fixed_cutoff.
	# For now, lets require that both timepoints are detected and not yet fixed.
	# There is an issue related to comparing fixed genotypes against non-fixed genotypes.
	left_was_fixed = widgets.fixed(left_trajectory, process.fixed_limit)
	right_was_fixed = widgets.fixed(right_trajectory, process.fixed_limit)

	if left_was_fixed == right_was_fixed:
		left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, process.detection_limit, process.fixed_limit,
			inner = False)
	else:
		left_reduced, right_reduced = widgets.get_valid_points(left_trajectory, right_trajectory, process.detection_limit, inner = False)

	if left_reduced.empty or right_reduced.empty:
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
	left_fixed: pandas.Series = left[left.gt(fixed_cutoff)].dropna()
	right_fixed: pandas.Series = right[right.gt(fixed_cutoff)].dropna()
	if left_fixed.empty or right_fixed.empty:
		# Both are undetected, since they both failed the invalid timepoint filter.
		p_value = False
	else:
		# Since these genotypes fixed immediately, they should only be grouped together if they
		# fixed at the same timepoint.
		fixed_at_same_time = left_fixed.index[0] == right_fixed.index[0]
		# Check if the trajectories overlap completely
		overlap = set(left_fixed.index) & set(right_fixed.index)
		complete_overlap = len(overlap) == len(left_fixed) and len(overlap) == len(right_fixed)
		p_value = fixed_at_same_time and complete_overlap
	result = 0 if p_value else math.nan

	return result


def plot_benchmark_results(benchmarks, output_filename):
	import matplotlib.pyplot as plt
	import seaborn
	labels = list()
	values = list()
	for k, v in benchmarks.items():
		labels.append(k)
		values.append(v)
	fig, ax = plt.subplots(figsize = (12, 10))
	seaborn.barplot(x = values, y = labels)

	plt.xlabel('time in seconds', fontsize = 14)
	plt.ylabel('number of processes', fontsize = 14)
	plt.title('Serial vs. Multiprocessing via Parzen-window estimation', fontsize = 18)
	plt.grid()

	plt.savefig(output_filename)


def benchmark(filename):
	import time

	table = pandas.read_excel(filename)
	table = table.set_index('Trajectory')
	table = table[widgets.get_numeric_columns(table.columns)]
	combinations = list(itertools.combinations(table.index, 2))
	benchmarks = dict()
	for index in [None] + list(range(1, 9)):
		logger.info(f"Running iteration {index}")
		distance_calculator = DistanceCalculator(.03, .97, 'binomial', threads = index)
		start = time.time()
		if index is None:
			results = [calculate_distance(distance_calculator, e, table) for e in combinations]
		else:
			pool = multiprocessing.Pool(processes = index)
			results = [pool.apply_async(calculate_distance, args = (distance_calculator, e, table)) for e in combinations]
			results = [i.get() for i in results]
		duration = time.time() - start
		label = str(index)
		benchmarks[label] = duration
	return benchmarks


if __name__ == "__main__":
	b = benchmark()
	plot_benchmark_results(b)
