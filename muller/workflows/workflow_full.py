#!python
"""
	Main script to run the muller workflow.
	.
	|---- .clade.annotated.(png|svg)
	|---- .genotypes.tsv
	|---- .lineage.clade.(png|svg)
	|---- figures/
	|----|---- clade/
	|----|----|---- .clade.annotated.(png|svg)
	|----|----|---- .clade.unannotated.(png|svg)
	|----|----|---- .clade.genotypes.(png|svg)
	|----|----|---- .clade.trajectories.(png|svg)
	|----|---- unique/
	|----|----|---- .unique.annotated.(png|svg)
	|----|----|---- .unique.unannotated.(png|svg)
	|----|----|---- .unique.genotypes.(png|svg)
	|----|----|---- .unique.trajectories.(png|svg)
	|----|---- .dendrogram.png
	|----|---- .pairwisedistance.svg
	|---- scripts/
	|----|---- .r
	|---- suplementary-files/
	|----|---- .genotypeinformation.json
	|----|---- .lineagescores.tsv
	|----|---- .options.json
	|---- tables/
	|----|---- .trajectories.tsv
	|----|---- .genotypes.tsv
	|----|---- .distance.tsv
	|----|---- .linkagematrix
	|----|---- .ggmuller.edges.tsv
	|----|---- .ggmuller.populations.tsv
	|----|---- .mullerdataframe.tsv
"""
import argparse
from pathlib import Path
from typing import List, Optional, Union

import pandas
from loguru import logger

from muller import clustering, dataio, inheritance
from muller.dataio import projectdata


def run_genotype_inference_workflow(trajectoryio: Union[str, Path, pandas.DataFrame], metric: str, dlimit: float, slimit:float, flimit: float, pvalue: float,
		breakpoints: List[float], known_genotypes: Optional[Path], threads: Optional[int]) -> projectdata.DataGenotypeInference:
	"""
	Parameters
	----------
	trajectoryio: Union[str,Path,pandas.DataFrame]
		Either the path to a table with mutational trajectories or the table itself.
		The table should already be cleaned/filtered.
		The table should be organized with timepoints as columns and rows as genotype frequencies.
		Ex.
		genotype	0	7	13 ...
		genotype-1	0	.1	.3 ...
		genotype-2	0	.2	.1 ...
	metric: Literal['']
	dlimit, slimit, flimit, pvalue: float
	breakpoints
	known_genotypes
	threads
	"""
	if isinstance(trajectoryio, (str, Path)):
		logger.info(f"Reading '{trajectoryio}' as the trajectory table.")
		trajectories, trajectory_info = dataio.parse_trajectory_table(trajectoryio)
	else:
		trajectories = trajectoryio
		trajectory_info = None

	# Read in the starting genotypes.
	if known_genotypes:
		known_genotypes = dataio.parse_known_genotypes(known_genotypes)

	genotype_generator = clustering.ClusterMutations(
		metric = metric,
		dlimit = dlimit,
		slimit = slimit,
		flimit = flimit,
		pvalue = pvalue,
		starting_genotypes = known_genotypes,
		threads = threads
	)

	genotype_data = genotype_generator.run(trajectories)

	return genotype_data


def run_genotype_lineage_workflow(genotypeio: Union[str, Path, pandas.DataFrame], dlimit: float, flimit: float, pvalue: float,
		known_ancestry: Optional[Path]) -> projectdata.DataGenotypeLineage:
	"""

	Parameters
	----------
	genotypeio
	dlimit
	flimit
	pvalue
	known_ancestry

	Returns
	-------

	"""
	lineage_generator = inheritance.LineageWorkflow(
		dlimit = dlimit,
		flimit = flimit,
		pvalue = pvalue
	)

	# Read in the input data if it is not already a pandas.DataFrame object
	if isinstance(genotypeio, (str, Path)):
		logger.info(f"Reading '{genotypeio}' as the genotype table.")
		genotypes = dataio.parse_genotype_table(genotypeio)
	else:
		genotypes = genotypeio

	lineage_data = lineage_generator.run(genotypes, known_ancestry)
	return lineage_data


def run_ggmuller():
	pass


def run_workflow(program_options: argparse.Namespace):
	trajectoryio = program_options.trajectories
	result_genotype_inference = run_genotype_inference_workflow(
		trajectoryio,
		program_options.metric,
		program_options.dlimit,
		program_options.flimit,
		program_options.pvalue,
		program_options.detection_breakpoint,
		program_options.breakpoints,
		program_options.starting_genotypes,
		program_options.threads
	)

	result_genotype_lineage = run_genotype_lineage_workflow(
		result_genotype_inference.table_genotypes, # should already be sorted.
		dlimit = program_options.dlimit,
		flimit = program_options.flimit,
		pvalue = program_options.pvalue,
		known_ancestry = program_options.known_ancestry
	)

if __name__ == "__main__":
	filename = Path("/home/cld100/Documents/github/muller_diagrams/tests/data/tables_input_trajectories/B1_muller_try1.trajectories.original.tsv")
	result = run_genotype_inference_workflow(
		trajectoryio = filename,
		metric = "binomialp",
		dlimit = 0.03,
		slimit = 0.15,
		flimit = 0.97,
		pvalue = 0.05,
		breakpoints = [],
		known_genotypes = None,
		threads = 8
	)
	print(result)
