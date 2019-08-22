from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger

try:
	from muller.inheritance import scoring
	from muller.inheritance.genotype_ancestry import Ancestry
	from muller import widgets
except ModuleNotFoundError:
	from . import scoring
	from .genotype_ancestry import Ancestry
	from .. import widgets


class LineageWorkflow:
	"""
		Orders genotypes by which background they belong to.
	Parameters
	----------
	dlimit: float
		The detection limit
	flimit: float
		The cutoff value to consider a genotype "fixed"
	additive_cutoff: float
		Used when testing whther a nested genotype is consistently greater than an unnested genotype
	subtractive_cutoff: float
		Used to test whether the combined frequencies of the nested/unnested genotype are consistently greater than the fixed cutoff.
	derivative_cutoff: float
		Used when testing whether two genotypes are correlated, not correlated, or anticorrelated. correlated/anticorrelated genotypes
		must have a covariance outside the range [-`derivative_cutoff`, `derivative_cutoff`].
	"""

	def __init__(self, dlimit: float, flimit: float, additive_cutoff: float, subtractive_cutoff: float, derivative_cutoff: float, ):
		self.dlimit = dlimit
		self.flimit = flimit
		self.additive_cutoff = additive_cutoff
		self.subtractive_cutoff = subtractive_cutoff
		self.derivative_cutoff = derivative_cutoff
		self.genotype_nests: Optional[Ancestry] = None

		#TODO: Add configurable weights for the score tests.

		self.score_records:List[Dict[str,float]] = list() # Keeps track of the individual score values for each pair

	def add_known_lineages(self, known_ancestry: Dict[str, str]):
		for identity, parent in known_ancestry.items():
			# TODO need to add a way to prevent circular ancestry links when a user manually assigns ancestry. Current workaround forces the manual parent to be in the root background.
			logger.debug(f"Adding {parent} as a potential background for {identity}")
			self.genotype_nests.add_genotype_to_background(parent, 'genotype-0', priority = 100)
			# Use a dummy priority so that it is selected before other backgrounds.
			self.genotype_nests.add_genotype_to_background(identity, parent, priority = 100)

	def score_pair(self, nested_genotype: pandas.Series, unnested_trajectory) -> Dict[str,float]:
		detected_left, detected_right = widgets.get_valid_points(nested_genotype, unnested_trajectory, dlimit = self.dlimit, inner = True)
		score_greater = scoring.calculate_greater_score(nested_genotype, unnested_trajectory, self.additive_cutoff)
		#score_subtractive = scoring.calculate_summation_score(nested_genotype, unnested_trajectory, self.flimit, self.subtractive_cutoff)
		score_subtractive = scoring.calculate_summation_score(detected_left, detected_right, self.flimit, self.subtractive_cutoff)
		score_area = scoring.calculate_area_score(nested_genotype, unnested_trajectory)

		total_score = score_greater + score_subtractive + score_area
		logger.debug(f"{nested_genotype.name}\t{unnested_trajectory.name}\t{score_greater}\t{score_subtractive}\t{score_area}")
		if total_score > 0:
			# The derivative check is only useful when deciding between possible candidates, since it does not provide evidence itself that a
			# genotype is a potential background. So, at least one of the other checks should have been passed with no
			# evidence against the candidate background.

			# The derivative score should only be computed using the timepoints where the series overlap.

			score_derivative = scoring.calculate_derivative_score(detected_left, detected_right, detection_cutoff = self.dlimit,
				cutoff = self.derivative_cutoff)
			# Note that a previous version accidentlly added the derivative cutoff to the total score.
			total_score += score_derivative
		else:
			score_derivative = None
		score_data = {
			'nestedGenotype': nested_genotype.name,
			'unnestedGenotype': unnested_trajectory.name,
			'scoreAdditive': score_greater,
			'scoreSubtractive': score_subtractive,
			'scoreArea': score_area,
			'scoreDerivative': score_derivative,
			'totalScore': total_score

		}
		return score_data

	def show_ancestry(self, sorted_genotypes: pandas.DataFrame):
		logger.log("COMPLETE", "Final Ancestry:")
		for genotype_label in sorted_genotypes.index:
			candidate = self.genotype_nests.get_highest_priority(genotype_label)
			logger.log('COMPLETE', f"{genotype_label}\t{candidate}")

	def run(self, sorted_genotypes: pandas.DataFrame, known_ancestry: Dict[str, str] = None) -> Ancestry:
		"""
			Infers the lineage from the given genotype table.
		Parameters
		----------
		sorted_genotypes: pandas.DataFrame
			A dataframe of sorted genotypes based on when the genotype was first detected and first fixed.
		known_ancestry: Dict[str,str]
			Manually-assigned ancestry values. For now, the parent genotype is automatically assigned to the root genotype to prevent
			circular links from forming.
		"""
		initial_background = sorted_genotypes.iloc[0]
		self.genotype_nests = Ancestry(initial_background, timepoints = sorted_genotypes)
		self.add_known_lineages(known_ancestry if known_ancestry else dict())

		for unnested_label, unnested_trajectory in sorted_genotypes[1:].iterrows():
			logger.debug(f"Nesting {unnested_label}")
			# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
			test_table = sorted_genotypes[:unnested_label].iloc[::-1]
			for nested_label, nested_genotype in test_table.iterrows():
				if nested_label == unnested_label: continue
				score_data = self.score_pair(nested_genotype, unnested_trajectory)
				self.score_records.append(score_data)
				self.genotype_nests.add_genotype_to_background(unnested_label, nested_label, score_data['totalScore'])

		self.show_ancestry(sorted_genotypes)
		return self.genotype_nests


def get_maximum_genotype_delta(genotype_deltas: List[Tuple[str, float]]) -> Tuple[str, float]:
	if genotype_deltas:
		correlated_label, correlated_delta = max(genotype_deltas, key = lambda s: s[1])
	else:
		correlated_label = "N/A"  # Shouldn't end up being used.
		correlated_delta = 0
	return correlated_label, correlated_delta
