from pathlib import Path
import pandas
from muller.clustering import filters # Should probably move this out of the `clustering` folder.
from typing import Optional

class FilterWorkflow:
	def __init__(self, dlimit:float, flimit:float, filter_consistency: float = 0.01, filter_single:bool = True, filter_startfixed:bool = True):

		# Save the general application variables.
		self.dlimit = dlimit # The detection cutoff value. Values below this threshold are considered undetected.
		self.flimit = flimit # The value by which a series is considered 'fixed'
		self.filter_consistency = filter_consistency # The amount by which a series should fluctuate to not be discarded.
		self.filter_single = filter_single # Whether to filter out genotypes consisteng of a single measured value.
		self.filter_startsfixed = filter_startfixed # Whether to fiter out series which are 'fixed' at the initial timepoint.

		# Set up the two filtering processes.
		self.filter_trajectory = filters.TrajectoryFilter(
			detection_cutoff = self.dlimit,
			fixed_cutoff = self.flimit,
			filter_consistency = self.filter_consistency,
			filter_single = self.filter_single,
			filter_startfixed = self.filter_startsfixed
		)
		self.filter_genotypes = filters.GenotypeFilter(
			detection_cutoff = dlimit, fixed_cutoff = flimit, frequencies = []
		)
	def run(self, table:pandas.DataFrame, output_filename:Optional[Path]):
		pass

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"--filename",
		help = "The table of mutationsal trajectories to infer genotypes from.",
		type = Path,
	)