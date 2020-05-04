import math
import statistics
from typing import Dict, List, Tuple

import pandas
import scipy.stats as stats
from loguru import logger

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

	def calculate_score_greater(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
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
		"""
		combined = left + right
		above_fixed = combined - self.flimit

		# Finally, check whether the frequency of this series is consistently above the `cutoff` value.
		# Timepoints that do not sum to greater than the fixed breakpoint will be negative.
		result = above_fixed.mean() > self.pvalue

		return int(result)


class Score:
	""" Refactored as a class so that all the dlimit,flimit,pvalue,etc variables don't have to be passed around"""

	def __init__(self, dlimit: float, flimit: float, pvalue: float, weights:Tuple[int,int,int,int] = (1, 2, 1, 2)):
		self.pvalue = pvalue
		self.dlimit = dlimit
		self.flimit = flimit
		self.slimit = 0.15

		# Keep this around as a fallback
		# Shapely has been having issues, so may need to fallback to the legacy area score.
		self.legacy_scorer = LegacyScore(pvalue, dlimit, flimit)
		self.weight_greater = weights[0]
		self.weight_above_fixed = weights[1]
		self.weight_derivative = weights[2]
		self.weight_jaccard = weights[3]

		self.debug = False

	def _get_growth_regions(self, series: pandas.Series) -> List[bool]:
		""" Returns True for regions of positive growth and fixed."""
		difference_series = series.diff()

		boolseries = [(i > self.flimit or j > 0) for i, j in zip(series.values, difference_series.values)]
		return boolseries

	def calculate_score_greater(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
		use_advanced = False
		series_overlap = widgets.overlap(nested_genotype, unnested_genotype, self.dlimit)
		if series_overlap == 0:
			return self.weight_greater * -1

		# THe t-test has to be corrected for the case where the two series do not completely overlap
		nested_genotype, unnested_genotype = widgets.get_valid_points(nested_genotype, unnested_genotype, self.dlimit)
		if use_advanced:
			raise NotImplementedError
		else:
			score = self.calculate_score_greater_basic(nested_genotype, unnested_genotype)

		return float(score)  # Cast to float so the dtypes are consistent

	def calculate_score_greater_basic(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> int:
		difference = nested_genotype - unnested_genotype

		nested_is_above_unnested = difference > self.dlimit
		unnested_is_above_nested = difference < -self.dlimit

		# Select the timepoints that correspond to whichever was greater.
		nested = difference[nested_is_above_unnested]
		unnested = difference[unnested_is_above_nested]

		portion_of_series = 0.5
		# The `cutoff` value is just half the total number of timepoints.
		cutoff = len(difference) * portion_of_series

		nested_timepoints = len(nested) # Number of timepoints where the nested genotype was greater
		unnested_timepoints = len(unnested) # Number of timepoints where the unnested genotype was greater.
		# This is adding up the difference in values. Basically the area under the curve.
		# Since each series already omits the timepoints not in that category the abs() method just
		# makes sure the differences are positive.
		nested_sum = abs(nested.sum())
		unnested_sum = abs(unnested.sum())

		# This is basically the mean of the difference between series.
		nested_mean = abs(nested.mean())
		unnested_mean = abs(unnested.mean())

		# Make sure that even if nested is generally larger than unnested, that unnested was not also greater than nested significantly.

		is_consistently_greater = (nested_timepoints > cutoff) and (nested_mean > self.dlimit)
		is_significantly_greater = nested_sum > self.slimit
		left = is_significantly_greater or is_consistently_greater

		is_consistently_less = (unnested_timepoints > cutoff) and (unnested_mean > self.dlimit)
		is_significantly_less = unnested_sum > self.slimit
		right = is_significantly_less or is_consistently_less

		if left and right:
			score = math.nan
		elif is_consistently_greater or is_significantly_greater:
			score = 1
		elif is_consistently_less or is_significantly_less:
			score = -1
		else:
			score = math.nan

		return score

	def _single_sample_ttest(self, left:pandas.Series, right:pandas.Series):
		""" Implements a single sample ttest."""
		combined_series = (left + right) - (1+self.dlimit)
		logger.debug(combined_series.tolist())
		# Test if the result is greater than 0
		statistic, pvalue = stats.ttest_1samp(combined_series.values, 0)
		return statistic, pvalue
	def _multiple_sample_ttest(self, left:pandas.Series, right:pandas.Series):
		combined_series = (left + right).tolist()[1:]

		mean_c = statistics.mean(combined_series)
		var_c = self.dlimit
		mean_f = 1 + self.dlimit
		var_f = self.dlimit

		forward_statistic, forward_pvalue = stats.ttest_ind_from_stats(
			mean1 = mean_c,
			std1 = var_c ** 2,
			nobs1 = len(combined_series),
			mean2 = mean_f,
			std2 = var_f ** 2,
			nobs2 = len(combined_series)
		)
		return forward_statistic, forward_pvalue

	def calculate_score_above_fixed(self, left: pandas.Series, right: pandas.Series) -> int:
		"""
			Tests whether two genotypes consistently sum to a value greater than the fixed breakpoint. This suggests that one of the genotypes
			is in the background of the other, since otherwise the maximum combined frequency should, at most, be equal to the fixed cutoff value.
			Keep in mind that the variance is defined as the uncertainty in the measurements rather than computed using the given values.
		Parameters
		----------
		left, right: pandas.Series
		"""
		# Including points where one genotype was not detected will skew the results.
		left, right = widgets.get_valid_points(left, right, dlimit = self.dlimit, inner = True)
		combined_series = (left + right).tolist()[1:]

		if len(combined_series) == 0:
			result = 0
		elif len(combined_series) == 1:
			result = combined_series[0] > self.flimit
		else:
			forward_statistic, forward_pvalue = self._multiple_sample_ttest(left, right)
			#forward_statistic, forward_pvalue = self._single_sample_ttest(left, right)
			# Since we're using a two-sided test we need to convert it to a one-sided test.
			result = forward_pvalue / 2 < self.pvalue and forward_statistic > 0

		return int(result)

	def calculate_score_area(self, nested_genotype: pandas.Series, unnested_genotype: pandas.Series) -> float:
		"""
			Calculates a score based on the probability the unnested genotype is a subset of the nested_genotype.
			This take into account both the area of the nested genotype and the area of all other genotypes not in the nested genotype.
		Parameters
		----------
		nested_genotype, unnested_genotype: pandas.Series

		"""

		# If the nested genotype is not fixed, group the remaining frequencies into an `other` category.
		difference_series = nested_genotype - unnested_genotype
		if difference_series.mean() > 0:
			# noinspection PyTypeChecker
			other_genotypes: pandas.Series = self.flimit - nested_genotype
		else:
			# noinspection PyTypeChecker
			other_genotypes: pandas.Series = self.flimit - unnested_genotype  # In case we're testing if a small genotype contains a large genotype

		other_genotypes = other_genotypes.mask(lambda s: s < 0, 0.0001)  # Since the flimit is not exactly 1.

		unnested_polygon = polygon.as_polygon(unnested_genotype)

		nested_polygon = polygon.as_polygon(nested_genotype)
		other_polygon = polygon.as_polygon(other_genotypes)

		is_subset_nested = areascore.is_subset_polygon(nested_polygon, unnested_polygon)
		is_subset_other = areascore.is_subset_polygon(other_polygon, unnested_polygon)
		is_subset_nested_reversed = areascore.is_subset_polygon(unnested_polygon, nested_polygon)  # Check the reverse case

		nested_area = areascore.area_of_series(nested_genotype)
		unnested_area = areascore.area_of_series(unnested_genotype)
		common_area_nested = areascore.X_and_Y_polygon(unnested_polygon, nested_polygon)
		xor_area_unnested = areascore.difference_polygon(unnested_polygon, nested_polygon)  # This does not distinguish between xor left vs xor right


		if self.debug:
			logger.debug(
				f"calculate_score_jaccard()->({is_subset_nested}, {is_subset_other}, {is_subset_nested_reversed}), ({common_area_nested:.2f}, {xor_area_unnested:.2f})")
		if is_subset_nested and is_subset_other:
			# Evidence for both scenarios
			# Test if the nested genotype is sufficiently large to assume the unnested genotype is a subset.
			# Test only the area where the unnested genotype was detected.
			other_polygon = polygon.as_polygon(other_genotypes)
			common_area_other = areascore.X_and_Y_polygon(unnested_polygon, other_polygon)
			score = int(common_area_nested > 2 * common_area_other)

		elif is_subset_nested:
			score = 1
		elif is_subset_nested_reversed and not is_subset_other:
			score = -1
		else:
			score = 0
		if score == 0 and xor_area_unnested > common_area_nested * 2:
			score = -1
		elif unnested_area > 2 * nested_area:
			score = -1
		score = score * self.weight_jaccard
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
		valid_left, valid_right = widgets.get_valid_points(left, right, self.dlimit, self.flimit, inner = True)

		if valid_left.empty:
			score = 0
		elif len(valid_left) > 20 or True:
			dotproduct, correlated_timepoints = self.derivative(valid_left, valid_right)
			#logger.debug(f"{dotproduct}, {correlated_timepoints}, {len(valid_left)}")
			#logger.debug(valid_left.diff().tolist())
			#logger.debug(valid_right.diff().tolist())
			if dotproduct > 0.01:
				score = 1
			elif dotproduct < -0.01:
				score = -1
			else:
				score = 0
		else:
			covariance = valid_left.cov(valid_right)
			if covariance > 0.01: score = 1
			elif covariance < -0.01: score = -1
			else: score = 0
		score = score * self.weight_derivative
		return score

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
			score_fixed = self.legacy_scorer.calculate_summation_score(detected_left, detected_right)
		else:
			score_fixed = self.calculate_score_above_fixed(detected_left, detected_right)

		score_greater = self.calculate_score_greater(detected_left, detected_right)
		if math.isnan(score_greater): score_greater = 0
		score_area = self.calculate_score_area(nested_genotype, unnested_trajectory)

		total_score = score_fixed + score_greater + score_area
		if self.debug:
			logger.debug(f"{nested_genotype.name}\t{unnested_trajectory.name}\t{score_fixed}\t{score_greater}\t{score_area}\t{total_score}")
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
			score_derivative = math.nan
		if math.isnan(score_derivative): score_derivative = 0
		score_data = {
			'nestedGenotype':   nested_genotype.name,
			'unnestedGenotype': unnested_trajectory.name,
			'scoreGreater':     score_greater,
			'scoreFixed':       score_fixed,
			'scoreArea':        score_area,
			'scoreDerivative':  score_derivative,
			'totalScore':       total_score
		}
		return score_data

	def derivative(self, left: pandas.Series, right: pandas.Series) -> Tuple[float, int]:
		# l, r = widgets.get_valid_points(left, right, 0.03, 0.97, inner = True)
		l, r = left, right
		normalize = lambda s: 1 if s > self.dlimit else (-1 if s < -self.dlimit else 0)

		ldiff = l.diff()[1:]
		rdiff = r.diff()[1:]

		result = ldiff.dot(rdiff)

		ldiff = ldiff.apply(normalize)
		rdiff = rdiff.apply(normalize)

		nresult = ldiff.dot(rdiff)
		return result, nresult
