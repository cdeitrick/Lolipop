from pathlib import Path
import pandas
from muller import dataio
from typing import Union
from loguru import logger
from muller import clustering
import argparse
from typing import List, Optional
class GenotypeInferenceWorkflow:
	""" Infers genotypes from the available trajectory frequencies
		This workflow expects that the data has already been cleaned and processed with desired filters.
	"""

	def __init__(self, program_options):
		logger.critical(f"Have to provide the frequency breakpoints.")
		breakpoints = []
		self.program_options = program_options
		self.genotype_generator = clustering.ClusterMutations(
			method = self.program_options.method,
			metric = self.program_options.metric,
			dlimit = self.program_options.detection_breakpoint,
			flimit = self.program_options.fixed_breakpoint,
			pvalue = self.program_options.pvalue,
			dbreakpoint = self.program_options.detection_breakpoint,
			breakpoints = breakpoints,
			starting_genotypes = self.program_options.starting_genotypes,
			trajectory_filter = self.trajectory_filter,
			genotype_filter = self.filter_genotype,
			filename_pairwise = self.program_options.filename_pairwise,
			threads = self.program_options.threads
		)

	def run(self, trajectoryio:Union[str,Path, pandas.DataFrame]):
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
		"""
		if isinstance(trajectoryio, (str, Path)):
			logger.info(f"Reading '{trajectoryio}' as the trajectory table.")
			trajectories = dataio.parse_genotype_table(trajectoryio)
		else:
			trajectories = trajectoryio

def create_parser(args:Optional[List] = None)->argparse.Namespace:
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"filename",
		help = "The timeseries table of the observed mutational trajectories.",
		type = Path,
		required = True
	)

if __name__ == "__main__":
	test_table = Path()