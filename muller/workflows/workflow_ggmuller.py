from pathlib import Path
from typing import Dict, Optional, Union

import pandas
from loguru import logger

from muller import dataio, graphics


# noinspection PyUnresolvedReferences
class GGMuller:
	""" Uses the `population` and `edges` table to generate muller plots and timeseries.
		Parameters
		----------
		populations: Union[Path. pandas.DataFrame]
		The population table is essentially a timeseries table in longformat (one observation per row).
		* `Generation`: The timepoint for this row.
		* `Identity`: Name of the series to which this observation belongs.
		* `Population` The measured frequency (as a number in [0,100])
		edges: Union[Path, pandas.DataFrame]
		Describes the lineage of the genotypes.
		* `Identity`: The genotype id.
		* `Parent`: The parent of `Identity`

		Output Files
		------------
		`data/population.tsv': The input population table
		`data/edges.tsv`: The input edges table
		`data/timeseries.tsv`: The condensed vertion of `population.tsv` where
			each column corresponds to a specific timepoint and each row corresponds to a
			single genotype.
		`figures/muller.(png|svg)` The muller diagram.

		Output Structure
		----------------
		folder
		|---- 'data'
		|----|---- 'population.tsv'
		|----|---- `edges.tsv'
		|----|---- `timeseries'
		|----'figures'
		|----|---- 'muller.(png|svg)

	"""

	def __init__(self, render: bool = True, folder: Optional[Path] = None):
		self.render = render
		if folder is None:
			self.paths = dict()
		else:
			self.paths = self.generate_filesystem_structure(folder)

	##############################################################################################################################################
	# --------------------------------------------------------------- Utilities ---------------------------------------------------------------------
	##############################################################################################################################################

	@staticmethod
	def generate_filesystem_structure(parent_folder: Path) -> Dict[str, Path]:
		parent_folder = dataio.checkdir(parent_folder)
		folder_data = dataio.checkdir(parent_folder / "data")
		folder_figure = dataio.checkdir(parent_folder / "figures")
		paths = {
			'outputTablePopulationLong': folder_data / "populations.long.tsv",
			'outputTablePopulationWide': folder_data / "populations.wide.tsv",
			'outputTableEdges':          folder_data / "edges.tsv",
			'outputTableMuller':         folder_data / "mullerformat.tsv",
			'outputFigureMuller':        folder_figure / "mullerplot.png"
		}
		return paths

	@staticmethod
	def load_table(io: Union[str, Path, pandas.DataFrame]) -> pandas.DataFrame:
		if not isinstance(io, pandas.DataFrame):
			# Assume it is one of the other formats
			result = dataio.import_table(io)
		else:
			result = io
		return result

	@staticmethod
	def convert_to_timeseries(table: pandas.DataFrame):
		result = table.pivot(index = 'Identity', columns = 'Generation', values = 'Population')
		return result

	##############################################################################################################################################
	# ------------------------------------------------------------------- main ---------------------------------------------------------------------
	##############################################################################################################################################

	def run(self, populations: Union[str, Path, pandas.DataFrame], edges: Union[str, Path, pandas.DataFrame], output_folder: Optional[Path] = None) -> \
	Dict[str, Path]:
		"""
		Parameters
		----------
		populations: Union[Path. pandas.DataFrame]
		edges: Union[Path, pandas.DataFrame]
		output_folder: Path
			The folder to put the output in. This will override the defaults generated when this workflow was initialized.

		Returns
		-------
		paths: Dict[str,Path]
			The dictionary with the predetermined filenames.
		"""
		if output_folder is not None:
			self.paths = self.generate_filesystem_structure(output_folder)
		# Load the two tables if they ar not already loaded.
		table_populations = self.load_table(populations)
		# The edges table maps genotypes to their corresponding parent.
		# Since this is a simple mapping, convert it to a pandas.Series object
		# 	which makes handling key-value pairs a little more convienient.
		table_edges = self.load_table(edges)
		# Need to save the pandas.Series object as a separate variable so that it isn't changed below.
		series_edges = table_edges.set_index('Identity')['Parent']
		# convert `table_populations` into a wide-format (columns == timepoints, rows = samples) timeseries table.
		table_populations_wide = self.convert_to_timeseries(table_populations)

		# Convert the `populations` and `edges` tables to a muller-df format.
		# The muller_df combines the information from those two tables.
		# The resulting dataframe will have the following columns:
		# Generation	Identity	Population	Frequency	Group_id	Unique_id
		# Generation
		# Identity: The name of the series
		# Groupid: Since each series is plotted twice (once for the bottom portion and again for the top portion)
		# 	whenever there are child genotypes, Each of those two series needs a way of identifying them.
		#	The general format is `identity` + (a|b)
		# Unique_id: Similar to the `Group_id` above, but for each individual timepoint.
		#	ex. `identity` + (a|b) + '_' + `timepoint`

		muller_formatter = dataio.GenerateMullerDataFrame()
		table_muller = muller_formatter.run(series_edges, table_populations)

		# Save the tables before rendering the muller plot in case something goes wrong.
		self.save_table(table_edges, "outputTableEdges")
		self.save_table(table_populations, 'outputTablePopulationLong')
		self.save_table(table_populations_wide, 'outputTablePopulationWide')
		self.save_table(table_muller, 'outputTableMuller')

		filename_figure_muller = self.paths['outputFigureMuller']

		# Generate the muller plot.
		diagram_generator = graphics.MullerPlot(outlines = True, render = True)
		diagram_generator.plot(table_muller, filename_figure_muller)

		return self.paths

	##############################################################################################################################################
	# --------------------------------------------------------------- Output ---------------------------------------------------------------------
	##############################################################################################################################################

	def save_table(self, table, output: Union[str, Path]) -> Path:
		""" Save the current table.
			Parameters
			----------
			table:pandas.DataFrame
				The source table to save. Can be any of the tables used in this workflow, the exact type of which is specified by `output`.
			output: Union{str, Path]
				The path to save the current table as. providing a `Path` object will save the table as that filename.
				If a key from `self.paths` is provided, the filename will lookup the corresponding default from `self.paths,

			Returns
			-------
			filename: Path
				The filename the table was saved as.

		"""
		if isinstance(output, Path):
			# No need to change anything.
			filename = output
		elif isinstance(output, str) and output in self.paths:
			# Should be a key from the `self.paths` dictionary.
			filename = self.paths[output]
		else:
			# `output` is not a Path or a key from `self.keys`.
			valid_keys = sorted(self.paths.keys())
			message = f"Expected a filename object or a key from {valid_keys}. Got '{output}' instead."
			logger.error(message)
			raise ValueError(message)
		table.to_csv(filename, sep = "\t")
		return filename


if __name__ == "__main__":
	pass
