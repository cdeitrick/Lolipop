import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas


@dataclass
class DataWorkflowBasic:
	# Used to organize the output from the workflow.DistanceCache
	# Should only cover relevant data about the scripts as a whole.
	version: str  # The version of the scripts
	filename: Path  # The filename of the input dataset.
	program_options: Any # The command-line arguments passed to the scripts.
	def save(self, folder:Path):
		options = vars(self.program_options)
		options['filename'] = self.filename
		options['version'] = self.version

		file_name = folder / "options.json"
		# Need to conver `Path` to `str`
		options = {key: (str(value)) for key, value in options.items()}
		file_name.write_text(json.dumps(options, sort_keys = True, indent = 4))
@dataclass
class DataGenotypeInference:
	"""
		Contains data related to trajectories and genotypes.

	"""
	# Include the original trajectory table.
	table_trajectories: pandas.DataFrame
	table_trajectories_info: Optional[pandas.DataFrame]
	# A table with inferred genotypes sorted by maximum frequency and earliest timepoint.
	# Each genotype is the mean of consituent trajectories.
	table_genotypes: pandas.DataFrame
	# Map of genotypes to member trajectories.
	genotype_members: Dict[str,List[str]]

	# A pairwise distance matrix showing the calculated distance between an given pair of trajectories.
	# Type: 'muller.clustering.metrics.distance_cache.DistanceCache'
	matrix_distance: Optional[Any]
	# Generated from the hierarchal clustering step. Links trajectories based on the pairwise distance.
	clusterdata: Optional["DataHierarchalCluster"]

	def save(self, folder:Path, prefix:str):
		"""
			Saves the data generated while infering genotypes.
			.folder
			|---- {prefix}.genotypes.tsv
			|---- {prefix}.trajectories.tsv
			|---- {prefix}.distancematrix.tsv
			|---- {prefix}.linkagetable.tsv
		"""
		delimiter = '\t'
		suffix = 'tsv'
		filename_table_trajectory = folder / (prefix + f'.trajectories.{suffix}')
		filename_table_genotypes = folder / (prefix + f'.genotypes.{suffix}')
		filename_table_distance_matrix = folder / (prefix + f'.distancematrix.{suffix}')
		filename_table_linkage_matrix = folder / (prefix + f'.linkagetable.{suffix}')

		# Include the `genotype_members` table in the `genotypes` table.
		# Convert the `genotype_members` from Dict[str,List[str]] to Dict[str,str]
		members = {key:'|'.join(values) for key, values in self.genotype_members.items()}
		# Add the genotype members to the `genotypes` table.
		self.table_genotypes['members'] = members

		# Save the data
		self.table_trajectories.to_csv(filename_table_trajectory, sep = delimiter)
		self.table_genotypes.to_csv(filename_table_genotypes, sep = delimiter)
		if self.clusterdata is not None:
			self.clusterdata.table_linkage.to_csv(filename_table_linkage_matrix, sep = delimiter)
		if self.matrix_distance is not None:
			self.matrix_distance.squareform().to_csv(filename_table_distance_matrix, sep = delimiter)

		# Need to remove the `members` column from the genotype table so that the graphics workflow uses a purely numeric table
		self.table_genotypes.pop('members')

@dataclass
class DataHierarchalCluster:
	""" Output from the hierarchal cluster. """
	clusters: List[List[str]]
	table_linkage: Optional[pandas.DataFrame]
	distance_cutoff: float
	distance_quantile: float
	def to_dict(self)->Dict[str,Any]:
		data = {
			'clusters': self.clusters,
			'distanceCutoff': self.distance_cutoff,
			'distanceQuantile': self.distance_quantile
		}

		return data


@dataclass
class DataGGmuller:
	""" Holds data related to the current implementation of ggmuller."""
	table_population: pandas.DataFrame
	table_edges: pandas.DataFrame
	script_r: str

	def save(self, folder:Path, prefix: str):
		"""
			.
			|---- .ggmuller.population.tsv
			|---- .ggmuller.edges.tsv
			|---- .script.r
		"""
		delimiter = '\t'
		suffix = 'tsv'

		filename_table_populations = folder / (prefix + f'.ggmuller.population.{suffix}')
		filename_table_edges = folder / (prefix + f'.ggmuller.edges.{suffix}')
		filename_script_r = folder / (prefix + f'.script.r')

		self.table_population.to_csv(filename_table_populations, sep = delimiter)
		self.table_edges.to_csv(filename_table_edges, sep = delimiter)
		filename_script_r.write_text(self.script_r)

@dataclass
class DataGenotypeLineage:
	""" Holds variables generated during the lineage inference step."""
	# The resulting scores of each unnested genotype to cancidate nested genotypes.
	table_scores: pandas.DataFrame
	clusters: Any # muller.inheritance.genotype_ancestry.Ancestry
	table_edges: pandas.Series
	table_populations: pandas.DataFrame
	table_muller: pandas.DataFrame

	def save(self, folder:Path, prefix:str):
		delimiter = "\t"
		suffix = "tsv"

		filename_table_scores = folder / (prefix + f'.lineage.scores.{suffix}')
		filename_table_populations = folder / (prefix + f".lineage.populations.{suffix}")
		filename_table_edges = folder / (prefix + f".lineage.edges.{suffix}")
		filename_table_muller = folder / (prefix + f".lineage.muller.{suffix}")

		self.table_scores.to_csv(filename_table_scores, sep = delimiter)
		self.table_populations.to_csv(filename_table_populations, sep = delimiter)
		self.table_edges.to_csv(filename_table_edges, sep = delimiter, header = True) # Added header option to avoid the FutureWarning
		self.table_muller.to_csv(filename_table_muller, sep = delimiter, index = False)


if __name__ == "__main__":
	pass