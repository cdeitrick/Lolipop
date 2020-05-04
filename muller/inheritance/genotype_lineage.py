from typing import Dict, List, Optional, Tuple

import pandas
from loguru import logger

try:
	from muller.inheritance import scoring
	from muller.inheritance.genotype_ancestry import Ancestry
	from muller import widgets, dataio
except ModuleNotFoundError:
	from . import scoring
	from .genotype_ancestry import Ancestry
	from .. import widgets, dataio


class LineageWorkflow:
	"""
		Orders genotypes by which background they belong to.
	Parameters
	----------
	dlimit: float
		The detection limit
	flimit: float
		The cutoff value to consider a genotype "fixed"
	pvalue: float
		The pvalue to use for statistical tests.
	"""

	def __init__(self, dlimit: float, flimit: float, pvalue: float, weights = (1, 1, 2, 2), conservative:bool = False,debug: bool = False):
		self.dlimit = dlimit
		self.flimit = flimit
		self.pvalue = pvalue
		self.debug = debug
		self.genotype_nests: Optional[Ancestry] = None
		self.conservative = conservative
		self.scorer = scoring.Score(self.dlimit, self.flimit, self.pvalue, weights)

	def __repr__(self)->str:
		return f"LineageWorkflow(dlimit = {self.dlimit}, flimit = {self.flimit}, pvalue = {self.pvalue})"

	def add_known_lineages(self, known_ancestry: Dict[str, str]):
		for identity, parent in known_ancestry.items():
			# TODO need to add a way to prevent circular ancestry links when a user manually assigns ancestry. Current workaround forces the manual parent to be in the root background.
			logger.debug(f"Adding {parent} as a potential background for {identity}")
			self.genotype_nests.add_genotype_to_background(parent, 'genotype-0', priority = 100)
			# Use a dummy priority so that it is selected before other backgrounds.
			self.genotype_nests.add_genotype_to_background(identity, parent, priority = 100)

	def show_ancestry(self, sorted_genotypes: pandas.DataFrame):
		for genotype_label in sorted_genotypes.index:
			candidate = self.genotype_nests.get_highest_priority(genotype_label)
			if self.debug:
				logger.debug(f"{genotype_label}\t{candidate}")

	def run(self, sorted_genotypes: pandas.DataFrame, known_ancestry: Dict[str, str] = None) -> dataio.projectdata.DataGenotypeLineage:
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
		self.genotype_nests = Ancestry(initial_background, timepoints = sorted_genotypes, cautious = self.conservative)
		self.add_known_lineages(known_ancestry if known_ancestry else dict())

		score_records: List[Dict[str, float]] = list()  # Keeps track of the individual score values for each pair

		for unnested_label, unnested_trajectory in sorted_genotypes[1:].iterrows():
			# Iterate over the rest of the table in reverse order. Basically, we start with the newest nest and iterate until we find a nest that satisfies the filters.
			test_table = sorted_genotypes[:unnested_label].iloc[::-1]
			for nested_label, nested_genotype in test_table.iterrows():
				if nested_label == unnested_label: continue
				score_data = self.scorer.score_pair(nested_genotype, unnested_trajectory)
				score_records.append(score_data)
				self.genotype_nests.add_genotype_to_background(unnested_label, nested_label, score_data['totalScore'])

		self.show_ancestry(sorted_genotypes)


		# Need to generate the population and edges tables.
		table_edges = self.genotype_nests.as_ancestry_table()
		population_table_generator = dataio.GGMuller(cutoff_detection = self.dlimit, adjust_populations = True)
		table_populations = population_table_generator.generate_ggmuller_population_table(
			table_edges,
			sorted_genotypes
		)

		# Need to generate the muller table
		muller_table_generator = dataio.GenerateMullerDataFrame()
		table_muller = muller_table_generator.run(table_edges, table_populations)

		output_data = dataio.projectdata.DataGenotypeLineage(
			table_scores = pandas.DataFrame(score_records),
			clusters = self.genotype_nests, # Used to extract the `edges` table.
			table_edges = self.genotype_nests.as_ancestry_table(),
			table_populations = table_populations,
			table_muller = table_muller
		)

		return output_data



def get_maximum_genotype_delta(genotype_deltas: List[Tuple[str, float]]) -> Tuple[str, float]:
	if genotype_deltas:
		correlated_label, correlated_delta = max(genotype_deltas, key = lambda s: s[1])
	else:
		correlated_label = "N/A"  # Shouldn't end up being used.
		correlated_delta = 0
	return correlated_label, correlated_delta
