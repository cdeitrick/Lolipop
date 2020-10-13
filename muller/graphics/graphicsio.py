from pathlib import Path
from typing import Dict, List

import pandas

from muller import graphics
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
	def __init__(self, project_folder: Path, palettes: List[Palette], sample_basename: str, render: bool = True, draw_outline: bool = True):
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

	@staticmethod
	def generate_dendrogram(linkage_matrix, distance_matrix, filename: Path) -> Path:
		# Only need the distance matrix fpr the labels
		labels = distance_matrix.squareform().index
		graphics.plot_dendrogram(linkage_matrix, labels, filename)
		return filename
	@staticmethod
	def generate_heatmap(squareform: pandas.DataFrame, filename: Path) -> Path:

		graphics.plot_heatmap(squareform, filename)
		return filename
