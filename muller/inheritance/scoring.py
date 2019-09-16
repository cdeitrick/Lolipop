import math
from typing import Dict, List

import pandas
import scipy.stats as stats
from loguru import logger
from shapely import geometry

try:
	from muller import widgets
	from muller.inheritance import areascore
	from muller.inheritance import polygon
except ModuleNotFoundError:
	from . import areascore
	from . import polygon


class LegacyScore:
	def __init__(self, pvalue: float, dlimit: float, flimit: float):
		self.pvalue = pvalue
		self.dlimit = dlimit
		self.flimit = flimit

	def calculate_greater_score(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
		"""
			Tests whether the nested genotype is consistenty larger than  the nested genotype. The scoring is as follows:
			`nested_genotype` always larger than `unnested_genotype`: 1
			`nested_genotype` consistently larger than `unnested_genotype`: 0.5
			`nested_genotype` not consistently larger than `unnested_genotype`: 0
			`nested_genotype` sometimes smaller than `unnested_genotype`: -0.5
		Parameters
		----------
		nested_genotype, unnested_genotype: pandas.Series
		cutoff: The minimum distance between each timepoint between `nested_genotype` and `unnested_genotype`
		"""
		difference = nested_genotype - unnested_genotype

		score = int(difference.sum() > (self.pvalue * len(difference) / 2))

		# There may be a point where the unnested genotype is greater than the nested genotype.

		if any(difference < 0):
			score -= 1
		# logger.debug(f"{unnested_genotype.name}\t{nested_genotype.name}\t{score}")
		return score

	def calculate_derivative_score(self, left: pandas.Series, right: pandas.Series) -> float:
		"""
			Tests whther the two series are correlated or anticorrelated with each other. The scoring is as follows:
			correlated: 2
			uncorrelated: 0
			anticorrelated: -2
		Parameters
		----------
		left, right: pandas.Series
			The two series to test.
		detection_cutoff: float
			Used to determine which points to use for the comparison.
		cutoff: float
			The cutoff value to determine whther the series are actually correlated/anticorrelated. Any value within the range [-`cutoff`, `cutoff`]
			results in a score of 0.
		"""
		# Pandas implementation of the derivative check, since it basically just checks for covariance.
		valid_left, valid_right = widgets.get_valid_points(left, right, self.dlimit, self.flimit)
		if valid_left.empty:
			covariance = math.nan
		else:
			covariance = valid_left.cov(valid_right)

		if covariance > 0.01: score = 2
		elif covariance < -0.01: score = -2
		else: score = 0
		return score

	@staticmethod
	def calculate_area(nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
		"""
			Calculates a score based on the probability the unnested genotype is a subset of the nested_genotype.
			This take into account both the area of the nested genotype and the area of all other genotypes not in the nested genotype.
		Parameters
		----------
		nested_genotype, unnested_genotype: pandas.Series
		"""

		# If the nested genotype is not fixed, group the remaining frequencies into an `other` category.
		# noinspection PyTypeChecker
		other_genotypes: pandas.Series = 1 - nested_genotype
		area_nested = areascore.area_of_series(nested_genotype)
		area_unnested = areascore.area_of_series(unnested_genotype)
		area_other = areascore.area_of_series(other_genotypes)
		common_area_nested = areascore.calculate_common_area(nested_genotype, unnested_genotype)
		common_area_other = areascore.calculate_common_area(other_genotypes, unnested_genotype)

		is_subset_nested = areascore.is_subset_legacy(area_nested, area_unnested, common_area_nested)
		is_subset_other = areascore.is_subset_legacy(area_other, area_unnested, common_area_other)

		if is_subset_nested and is_subset_other:
			score = int(common_area_nested > 2 * common_area_other)
		elif is_subset_nested:
			score = 2
		else:
			score = -2
		# logger.debug(f"{unnested_genotype.name}\t{nested_genotype.name}\t{score}")
		return score

	def calculate_summation_score(self, left: pandas.Series, right: pandas.Series) -> int:
		"""
			Tests whether two genotypes consistently sum to a value greater than the fixed breakpoint. This suggests that one of the genotypes
			is in the background of the other, since otherwise the maximum combined frequency should, at most, be equal to the fixed cutoff value.
		Parameters
		----------
		left, right: pandas.Series
		fixed_cutoff: The fixed breakpoint.
		cutoff: Governs whether two series consistently sum to a value greater than `fixed_cutoff`. Should be in the range [0,1]
		"""
		# TODO: Maybe the cutoff should be the variance as calculated from the clustering step.
		combined = left + right
		above_fixed = combined - self.flimit

		# Finally, check whether the frequency of this series is consistently above the `cutoff` value.
		# Timepoints that do not sum to greater than the fixed breakpoint will be negative.
		result = above_fixed.mean() > self.pvalue

		return int(result)


class Score:
	""" Refactored as a class so that all the dlimit,flimit,pvalue,etc variables don't have to be passed around"""

	def __init__(self, dlimit: float, flimit: float, pvalue: float, debug:bool = False):
		self.pvalue = pvalue
		self.dlimit = dlimit
		self.flimit = flimit

		# Keep this around as a fallback
		# Shapely has been having issues, so may need to fallback to the legacy area score.
		self.legacy_scorer = LegacyScore(pvalue, dlimit, flimit)

		self.weight_greater = 1
		self.weight_above_fixed = 1
		self.weight_derivative = 2
		self.weight_jaccard = 2

		self.debug = debug

	def _get_growth_regions(self, series: pandas.Series) -> List[bool]:
		""" Returns True for regions of positive growth and fixed."""
		difference_series = series.diff()

		boolseries = [(i > self.flimit or j > 0) for i, j in zip(series.values, difference_series.values)]
		return boolseries

	def calculate_score_greater(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> int:
		"""
			Tests whether the nested genotype is consistenty larger than  the nested genotype. The scoring is as follows:
			`nested_genotype` always larger than `unnested_genotype`: 1
			`nested_genotype` consistently larger than `unnested_genotype`: 0.5
			`nested_genotype` not consistently larger than `unnested_genotype`: 0
			`nested_genotype` sometimes smaller than `unnested_genotype`: -0.5
		Parameters
		----------
		nested_genotype, unnested_genotype: pandas.Series
		"""

		series_overlap = widgets.overlap(nested_genotype, unnested_genotype, 0.03)
		if series_overlap == 0:
			return self.weight_greater * -1


		forward_statistic, forward_pvalue = stats.ttest_rel(nested_genotype.values, unnested_genotype.values)
		forward_result = forward_pvalue / 2 < self.pvalue and forward_statistic > 0

		# Test whether the unnested genotype is greater than the nested genotype.
		reverse_statistic, reverse_pvalue = stats.ttest_rel(unnested_genotype.values, nested_genotype.values)
		reverse_result = reverse_pvalue / 2 < self.pvalue and reverse_statistic > 0

		if forward_result and not reverse_result:
			score = 1
		elif not forward_result and reverse_result:
			score = -1
		else:
			score = 0

		return self.weight_greater * score

	def calculate_score_above_fixed(self, left: pandas.Series, right: pandas.Series) -> int:
		"""
			Tests whether two genotypes consistently sum to a value greater than the fixed breakpoint. This suggests that one of the genotypes
			is in the background of the other, since otherwise the maximum combined frequency should, at most, be equal to the fixed cutoff value.
			Keep in mind that the variance is defined as the uncertainty in the measurements rather than computed using the given values.
		Parameters
		----------
		left, right: pandas.Series
		"""

		combined_series = (left + right).tolist()[1:]
		mean_c = sum(combined_series) / len(combined_series)
		var_c = self.dlimit

		mean_f = self.flimit
		var_f = self.dlimit

		forward_statistic, forward_pvalue = stats.ttest_ind_from_stats(
			mean1 = mean_c,
			std1 = var_c ** 2,
			nobs1 = len(combined_series),
			mean2 = mean_f,
			std2 = var_f ** 2,
			nobs2 = len(combined_series))

		forward_result = forward_pvalue / 2 < self.pvalue and forward_statistic > 0

		return int(forward_result)

	def calculate_score_area(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
		"""
			Calculates a score based on the probability the unnested genotype is a subset of the nested_genotype.
			This take into account both the area of the nested genotype and the area of all other genotypes not in the nested genotype.
		Parameters
		----------
		nested_genotype, unnested_genotype: pandas.Series

		"""

		# If the nested genotype is not fixed, group the remaining frequencies into an `other` category.
		# noinspection PyTypeChecker

		# other_genotypes: pandas.Series = self.flimit - nested_genotype.apply(lambda s: nested_genotype.mean())

		other_genotypes: pandas.Series = self.flimit - nested_genotype
		other_genotypes = other_genotypes.mask(lambda s: s < 0, 0)  # Since the flimit is not exactly 1.

		unnested_polygon = polygon.as_polygon(unnested_genotype)
		logger.debug(unnested_genotype.index)

		nested_polygon = polygon.as_polygon(nested_genotype)
		other_polygon = polygon.as_polygon(other_genotypes)
		logger.debug(nested_polygon)
		logger.debug(unnested_polygon)
		is_subset_nested = areascore.is_subset_polygon(nested_polygon, unnested_polygon)
		is_subset_other = areascore.is_subset_polygon(other_polygon, unnested_polygon)
		is_subset_nested_reversed = areascore.is_subset_polygon(unnested_polygon, nested_polygon)  # Check the reverse case

		common_area_nested = areascore.X_and_Y_polygon(unnested_polygon, nested_polygon)
		xor_area_unnested = areascore.difference_polygon(unnested_polygon, nested_polygon)

		if is_subset_nested and is_subset_other:
			# Evidence for both scenarios
			# Test if the nested genotype is sufficiently large to assume the unnested genotype is a subset.
			# Test only the area where the unnested genotype was detected.
			other_polygon = polygon.as_polygon(other_genotypes)
			common_area_other = areascore.X_and_Y_polygon(unnested_polygon, other_polygon)
			score = int(common_area_nested > 2 * common_area_other)

		elif is_subset_nested:
			score = 2
		elif is_subset_nested_reversed:
			score = -2
		else:
			score = 0

		if common_area_nested < 2 * xor_area_unnested:
			score = -2

		return score

	def calculate_score_derivative(self, left: pandas.Series, right: pandas.Series) -> float:
		"""
			Tests whther the two series are correlated or anticorrelated with each other. The scoring is as follows:
			correlated: 2
			uncorrelated: 0
			anticorrelated: -2
		Parameters
		----------
		left, right: pandas.Series
			The two series to test.
		"""
		# Pandas implementation of the derivative check, since it basically just checks for covariance.
		valid_left, valid_right = widgets.get_valid_points(left, right, self.dlimit, self.flimit)

		if valid_left.empty:
			covariance = math.nan
		else:
			covariance = valid_left.cov(valid_right)
		# TODO: Need to change this since the covariance is messed up when a series is removed from the population (neg growth)


		if covariance > 0.01: score = 2
		elif covariance < -0.01: score = -2
		else: score = 0
		return score

	def calculate_score_derivative_legacy(self, left, right):
		"""
			Tests whther the two series are correlated or anticorrelated with each other. The scoring is as follows:
			correlated: 2
			uncorrelated: 0
			anticorrelated: -2
		Parameters
		----------
		left, right: pandas.Series
			The two series to test.
		"""
		# Pandas implementation of the derivative check, since it basically just checks for covariance.
		valid_left, valid_right = widgets.get_valid_points(left, right, self.dlimit, self.flimit, inner = True)

		if valid_left.empty:
			covariance = math.nan
		else:
			covariance = valid_left.cov(valid_right)
		# TODO: Need to change this since the covariance is messed up when a series is removed from the population (neg growth)
		derivative_left = valid_left.diff()
		derivative_right = valid_right.diff()
		total = 0
		for l, r in zip(derivative_left, derivative_right):
			both_positive = (l>0) and (r>0)
			both_negative = (l<0) and (r<0)
			total += int(both_positive or both_negative)
		if right.name == 'genotype-orchid':
			logger.warning(f"Derivative score:{left.name} {total}, {len(valid_left)}")
		if total > len(valid_left)*.67:
			return 2
		elif total < len(valid_left) * .33:
			return -2
		else:
			return 0

	def score_pair(self, nested_genotype: pandas.Series, unnested_trajectory) -> Dict[str, float]:
		detected_left, detected_right = widgets.get_valid_points(nested_genotype, unnested_trajectory, dlimit = self.dlimit, inner = False)
		if self.debug:
			logger.debug(f"Scoring a pair of series:")
			logger.debug(f"\t{nested_genotype.name}\t{nested_genotype.values}\t{nested_genotype.index}")
			logger.debug(f"\t{unnested_trajectory.name}\t{unnested_trajectory.values}\t{unnested_trajectory.index}")
			logger.debug(f"{self.dlimit}, {self.flimit}, {self.pvalue}")
			logger.debug(f"The detected portion of the series: ")
			logger.debug(f"\t{detected_left.values}")
			logger.debug(f"\t{detected_right.values}")

		if len(detected_left) < 3:
			score_greater = self.legacy_scorer.calculate_summation_score(detected_left, detected_right)
		else:
			score_greater = self.calculate_score_above_fixed(detected_left, detected_right)

		score_subtractive = self.calculate_score_greater(detected_left, detected_right)
		score_area = self.calculate_score_area(nested_genotype, unnested_trajectory)

		total_score = score_greater + score_subtractive + score_area
		if self.debug:
			logger.debug(f"{nested_genotype.name}\t{unnested_trajectory.name}\t{score_greater}\t{score_subtractive}\t{score_area}\t{total_score}")
		if total_score > 0:
			# The derivative check is only useful when deciding between possible candidates, since it does not provide evidence itself that a
			# genotype is a potential background. So, at least one of the other checks should have been passed with no
			# evidence against the candidate background.

			# The derivative score should only be computed using the timepoints where the series overlap.
			detected_left, detected_right = widgets.get_valid_points(nested_genotype, unnested_trajectory, dlimit = self.dlimit, inner = True)
			score_derivative = self.calculate_score_derivative(detected_left, detected_right)
			# Note that a previous version accidentlly added the derivative cutoff to the total score.
			total_score += score_derivative
		else:
			score_derivative = None
		score_data = {
			'nestedGenotype':   nested_genotype.name,
			'unnestedGenotype': unnested_trajectory.name,
			'scoreAdditive':    score_greater,
			'scoreSubtractive': score_subtractive,
			'scoreArea':        score_area,
			'scoreDerivative':  score_derivative,
			'totalScore':       total_score
		}
		return score_data
