from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas
from loguru import logger

from muller import dataio, graphics
from muller.graphics import Palette

# The filenames should be mapped to their corresponding palette.
# The dictionary maps palette names to the filenames that should be generated.
GraphicPath = Dict[str, List[Path]]


class OutputGeneratorGraphics:
	""" Implemented to avoid having a `god object` controlling all output. It also allows us to decouple the logic
		used to generate the graphics from the logic used to save the data tables.

		This class is designed to be used in parts rather than as a whole so that any workflow (suck as the `ggmuller` workflow)
		can use only the methods it needs.

		There are only two variables needed for multiple processes: the palettes and whether to generate renders of data as well.


	"""
	# TODO: Need to rafactor this so it doesn't need the sample basename.
	def __init__(self, project_folder: Path, palettes: List[Palette], sample_basename:str, render:bool = True, draw_outline:bool = True):
		"""
			Use this method to set up some variables used by multiple workflows, such as `render`.
		"""
		# All we need from the basic_data is the base filename, so that should be passed directly.
		self.render = render
		self.extensions = [".png", ".svg"]

		# Coerce the given palettes into a list of Palette objects. These palettes will be visible by
		# many of the implemented methods.
		if isinstance(palettes, Palette):
			# In case only one Palette object was passed and was not already in a list.
			self.palettes = [palettes]
		else:
			self.palettes: List[Palette] = palettes

		# Each of the files will be prefixed with the population name.
		# Use dataio.PathsGraphics to set up default filenames.

		# Set up the individual graphic generators.
		self.generator_panel_timeseries = graphics.TimeseriesPanel(render = self.render)
		self.generator_plot_timeseries = graphics.TimeseriesPlot(render = self.render)
		self.generator_plot_muller = graphics.MullerPlot(outlines = draw_outline, render = self.render)

	##############################################################################################################################################
	# ---------------------------------------------------------- General Utilites ----------------------------------------------------------------
	##############################################################################################################################################

	@staticmethod
	def get_base_filename(data_basic) -> str:
		if data_basic.program_options.name:
			base_filename = data_basic.program_options.name
		else:
			base_filename = data_basic.filename.stem
		if data_basic.program_options['sheetname'] and data_basic.program_options['sheetname'] != 'Sheet1':
			base_filename += '.' + data_basic.program_options['sheetname']

		return base_filename

	##############################################################################################################################################
	# --------------------------------------------------------------- Panels ---------------------------------------------------------------------
	##############################################################################################################################################

	def generate_timeseries_panels_multiple(self, table_genotypes: pandas.DataFrame, table_trajectories: pandas.DataFrame,
			filenames: Dict[str, List[Path]] = None) -> Dict[str, List[Path]]:
		"""
			Generates a panel consisting of the inferred genotypes with an inset of the corresponding trajectories.
		Parameters
		----------
		table_genotypes, table_trajectories: Dataframes of the genotypes and trajectories to plot.
		filenames: Dict[str,List[str]]
			An optional list of filenames to save the panels as. Should be a dictionary mapping a colorscheme to a list of filenames.
			The filenames should be a list of paths so that multiple render types can be included.

		Returns
		-------
			List[Path]: A list of files that were generated.
		"""
		logger.info("Generating series plots...")
		# The timeseries panels should include the filtered trajectories in the trajectories inset plot.
		filenames = self.parse_filenames(filenames, default = self.paths.files_panel)


		return filenames

	def generate_timeseries_panels(self, table_genotypes:pandas.DataFrame, table_trajectories:pandas.DataFrame, filename:Path):
		"""
			Generates a panel consisting of the inferred genotypes with an inset of the corresponding trajectories.
		Parameters
		----------
		table_genotypes, table_trajectories: Dataframes of the genotypes and trajectories to plot.
		"""

		self.generator_panel_timeseries.run(
			table_genotypes,
			filenames = filename,
			timeseries_trajectory = table_trajectories
		)

	def generate_timeseries_trajectories(self, trajectories: pandas.DataFrame, trajectory_palette: Dict[str, str],
			filenames: Dict[str, List[Path]] = None) -> Path:
		"""
			Generates a plot of the trajectory frequencies over time.
			Parameters
			----------
			trajectories: pandas.DatFrame
				A dataframe where each column represents a timepoint and each row is a unique trajectory,
			trajectory_palette: Data[str,str]
				The palette to apply when plotting the timeseries.
			filenames: List[Path]
				An optional list of filenames to save the plot as. Overrides `GraphicsPath.timeseries_plot_trajectory()`.

			Returns
			-------
			The filename of the resulting plot.
		"""
		# The inset of the the trajectory plots in the above plots can be a little hard to read.
		# Generate a plot with only the trajectories. Not sure if its worth including both palette versions.
		filenames = self.parse_filenames(filenames, default = self.paths.filename_timeseries_trajectory)
		if filenames is None:
			filenames = self.paths.filename_timeseries_trajectory

		for palette in self.palettes:
			fnames = filenames[palette.name]
			self.generator_plot_timeseries.plot(
				timeseries = trajectories,
				palette = trajectory_palette,
				filenames = fnames
			)
		return filenames

	def save_graphics(self, muller_df: pandas.DataFrame, genotypes: dataio.GenotypeCollection, squareform: pandas.DataFrame,
			filtered_trajectories: pandas.DataFrame):


		self.generate_heatmap(squareform)
		self.generate_muller_plots(muller_df, genotypes.get('annotations'))
		self.generate_timeseries_panels(genotypes, filtered_trajectories)
		self.generate_geneology_plots(genotypes)

	def generate_dendrogram(self, linkage_matrix, distance_matrix, filename: Path) -> Path:
		# Only need the distance matrix fpr the labels
		labels = distance_matrix.squareform().index
		graphics.plot_dendrogram(linkage_matrix, labels, filename)
		return filename

	def generate_heatmap(self, squareform: pandas.DataFrame, filename: Optional[Path] = None) -> Path:
		if filename is None:
			filename = self.paths.distance_heatmap
		graphics.plot_heatmap(squareform, filename)
		return filename

	def generate_geneology_plot(self, genotypes: dataio.GenotypeCollection, table_edges: pandas.DataFrame, filename: GraphicPath = None) -> None:
		""" Generate the lineage plots."""
		logger.info("Generating Lineage Plots...")
		filenames = self.parse_filenames(filenames, self.paths.files_lineage_image)
		genotype_annotations = genotypes.get('annotations')
		for palette in self.palettes:
			graphics.flowchart(table_edges, palette, annotations = genotype_annotations, filenames = filenames)

	def generate_muller_plots(self, muller_df: pandas.DataFrame, genotype_annotations, filename: Optional[GraphicPath] = None):
		"""
			Generates a muller plot from the 'muller_df` object.
		Parameters
		----------
		muller_df
		genotype_annotations
		filenames: Dict[str,List[str]]
			Maps colorschemes with filenames.

		Returns
		-------

		"""
		logger.info("Generating muller plots...")
		filenames = self.parse_filenames(filenames, self.paths.files_muller_diagram)
		# Draw the muller diagrams
		for palette in self.palettes:
			filenames = self.paths.files_muller_diagram[palette.name]
			self.generator_plot_muller.plot(muller_df, filenames, palette, genotype_annotations)

	def validate_filenames(self, filenames: Any):
		""" Makes sure that the each colorscheme is mapped to a list of filesnames."""
		palette_names = [p.name for p in self.palettes]
		if not isinstance(filenames, dict):
			# Check whether the filenames varaible is formatted correctly.
			message = f"The filenames should be a dictionary mapping a palette to a list of available filenames."
		elif set(filenames.keys()) != set(palette_names):
			message = f"The palettes specified in the `filenames` variable ({filenames.keys()}) do not correspond to those in the `palettes` variable ({palette_names})"
		else:
			message = None
		if message:
			raise ValueError(message)

	def parse_filenames(self, filenames: Any, default: GraphicPath) -> GraphicPath:
		""" The filenames should correspond to a specific palette."""
		if filenames is None:
			filenames = default
		elif isinstance(filenames, (str, Path)):
			filenames = {'unique': [filenames]}
		filenames = {key: [Path(v) for v in value] for key, value in filenames.items()}
		self.validate_filenames(filenames)
		return filenames
