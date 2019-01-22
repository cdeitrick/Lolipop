from typing import Any, List, Tuple

import pandas

try:
	from muller_genotypes.metrics.pairwise_calculation_cache import PairwiseCalculation
	from muller_genotypes.methods.hierarchical_method import hierarchical_method
	from muller_genotypes.methods.twostep_method import twostep_method
except ModuleNotFoundError:
	from ..metrics.pairwise_calculation_cache import PairwiseCalculation
	from .hierarchical_method import hierarchical_method
	from .twostep_method import twostep_method


def calculate_genotypes_from_given_method(timepoints: pandas.DataFrame, pairwise_calculations: PairwiseCalculation, method: str,
		similarity_breakpoint: float, difference_breakpoint: float, starting_genotypes: List[List[str]]) -> Tuple[List[List[str]], Any]:
	"""
		Calculates genotypes for the population using the method defined by `method`
	Parameters
	----------
	timepoints: pandas.DataFrame
		A dataframe of trajectory timepoints.
	pairwise_calculations: PairwiseCalculation
		Contains the relevant calculations/metrics required by each method.
	method: {'hierarchy', 'matlab'}
		Specifies the clustering method to use.
	similarity_breakpoint: float
	difference_breakpoint:float
	starting_genotypes:List[List[str

	Returns
	-------
	Tuple[pandas.DataFrame, Any]
	"""
	if method == "matlab":
		genotypes = twostep_method(timepoints, pairwise_calculations, similarity_breakpoint, difference_breakpoint, starting_genotypes)
		linkage_matrix = None
	elif method == "hierarchy":
		genotypes, linkage_matrix = hierarchical_method(pairwise_calculations, similarity_breakpoint)
	else:
		raise ValueError(f"Invalid clustering method: {method}")
	return genotypes, linkage_matrix
