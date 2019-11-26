from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple
from muller.clustering.metrics.distance_cache import DistanceCache

import pandas

@dataclass
class DataWorkflowBasic:
	# Used to organize the output from the workflow.
	# Should only cover relevant data about the scripts as a whole.
	version: str  # The version of the scripts
	filename: Path  # The filename of the input dataset.
	program_options: Any # The command-line arguments passed to the scripts.

@dataclass
class DataGenotypeInference:
	"""
		Contains data related to trajectories and genotypes.
	"""
	#info: Optional[pandas.DataFrame]
	# A table with inferred genotypes sorted by maximum frequency and earliest timepoint.
	# Each genotype is the mean of consituent trajectories.
	table_genotypes: pandas.DataFrame
	# Map of genotypes to member trajectories.
	genotype_members: Dict[str,List[str]]
	# This is saved since it contains both filtered and unfiltered trajectories.
	original_trajectories: Optional[pandas.DataFrame]

	# A pairwise distance matrix showing the calculated distance between an given pair of trajectories.
	distance_matrix: Optional[DistanceCache]
	# Generated from the hierarchal clustering step. Links trajectories based on the pairwise distance.
	linkage_matrix: Optional[pandas.DataFrame]

	# A list of trajectories that were filtered out.
	filter_cache: List = field(default_factory =list)

@dataclass
class DataGenotypeLineage:
	""" Holds variables generated during the lineage inference step."""
	# The resulting scores of each unnested genotype to cancidate nested genotypes.
	score_history: List[Dict[str, float]]

	clusters: Any # muller.inheritance.genotype_ancestry.Ancestry


if __name__ == "__main__":
	pass