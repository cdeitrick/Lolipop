from typing import Dict

import matplotlib.pyplot as plt

try:
	from muller.graphics import MullerPlot, TimeseriesPlot, palettes, Palette
except (ModuleNotFoundError, ImportError):
	from .muller_plot import MullerPlot
	from .plottimeseries import TimeseriesPlot
	from . import palettes
	from .palettes import Palette

plt.clf()
# plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Optional, Tuple, List

from loguru import logger


class TimeseriesPanel:
	""" Plots the genotype frequencies over time with an inset of the original mutaitonal trajectories.
		Parameters
		----------
			render: bool
				Whether to generate additional renders in svg format of the resulting genotype plot.
			legend: bool
				Whether to include a legend in the plot.
	"""

	def __init__(self, render: bool, legend: bool = True):
		self.render = render

		self.filetypes = ['.png']
		if self.render:
			self.filetypes += ['.svg']

		# Parameters concerning the plot legend.
		self.legend = legend
		self.legend_font_properties = {'size': 12}  # The size of the font used to label each series in the legend.
		self.legend_location = 'right'
		self.legend_title = 'Genotypes'

		self.dpi = 250

		# Set the scaling factors.
		# The plots seem to look well up to 20 or so series, but should probably be increased after that.
		# The current values seem to work well with datasets of 10 to 20 genotypes and 15 to 20 timepoints.
		self.scale_max = 3
		self.reference_x = 20
		self.reference_y = 20

	def _get_scaling_factor(self, width_x: int, width_y: int) -> Tuple[int, int]:
		""" Since very large datasets require the plots to scale up, try to infer a decent scaling factor
			So that the resulting plots are actually legible.
		"""
		scale_x = round(width_x / self.reference_x)
		scale_y = round(width_y / self.reference_y)
		# Make sure they aren't 0
		if scale_x == 0: scale_x = 1
		if scale_y == 0: scale_y = 1

		if scale_x > self.scale_max:
			logger.debug(f"Timeseries Plot: Reducing x scale from {scale_x} to {self.scale_max}.")
			scale_x = self.scale_max
		if scale_y > self.scale_max:
			logger.debug(f"Timeseries Plot: Reducing y scale from {scale_y} to {self.scale_max}.")
			scale_y = self.scale_max
		return scale_x, scale_y


	def run(self, timeseries_genotype: pandas.DataFrame, palette: Palette, filename: Path,
			timeseries_trajectory: Optional[pandas.DataFrame] = None):
		"""
			Plots the clustered muller_genotypes.
		Parameters
		----------
		timeseries_genotype: pandas.DataFrame
			A dataframe where columns correspond to timepoints and rows correspond to individual series to plot.
			This is for a generic timeseries as described and isn't meant for a specific table.
		palette: Palette
			Has information related to how the trajectories relate to each genotype.
		filename: Optional[str, Path]
			Path to save the genotype plot. May include multiple filenames.
		timeseries_trajectory: pandas.DataFrame
			The table with each trajectory used to generate the genotypes.
		"""
		palette_genotypes = palette
		palette_trajectories = palette_genotypes.get_trajectory_palette()

		# Need to scale the graph so very large datasets are still easily distinguishable.
		number_of_timepoints = len(timeseries_genotype.columns)
		number_of_series = len(timeseries_genotype)

		# Calculate what the scale should be given the dataset shape. Any dataset smaller than the reference is automatically converted to 1.
		scale_x, scale_y = self._get_scaling_factor(number_of_timepoints, number_of_series)
		# Create the plotting object

		# Set up the plotting area. The genotype plot will cover all three columns while the trajectory plot will be relegated to the
		# Upper left. The legend will occupy the remaining space in the upper right.
		plotter = TimeseriesPlot()
		# Manually define the figure to control the figure size
		figure = plt.Figure(figsize = (12 * scale_x, 10 * scale_y))
		grid = plt.GridSpec(2, 3, wspace = 0.4, hspace = 0.3, figure = figure)

		# Plot all trajectories
		taxes: plt.Axes = plt.subplot(grid[0, 0:2])
		plotter.plot(timeseries_trajectory, palette_trajectories, taxes)
		gaxes: plt.Axes = plt.subplot(grid[1, :])
		# Plot clustered genotypes. Should be same as above, but colored based on genotype cluster.
		# Plot the mean of each cluster
		plotter.plot(timeseries_genotype, palette_genotypes, gaxes)

		if self.legend:
			bbox = self.get_legend_size(number_of_series)
			if bbox:
				plt.legend(
					loc = self.legend_location,
					ncol = 2,
					bbox_to_anchor = bbox,
					prop = self.legend_font_properties,
					title = self.legend_title,
				)
		if filename:
			self.save_figure(filename)

	@staticmethod
	def get_legend_size(number_of_genotypes: int) -> Optional[Tuple[float, float]]:
		if number_of_genotypes < 20:
			bbox = (1, 2)
		elif number_of_genotypes < 30:
			bbox = (1, 1.8)
		else:
			bbox = None

		return bbox

	def save_figure(self, filename: Path):
		""" Saves the diagram in the format specified by the filename's suffix"""
		plt.savefig(filename, dpi = self.dpi)


class MullerPanel:
	""" A panel of graphics describing the lineage of a dataset."""

	def __init__(self, vertical: bool = True, render: bool = True):
		self.vertical = vertical  # The orientation of the plots.
		self.vertical = False
		self.render = render

	def plot(self, timeseries: pandas.DataFrame, muller_df: pandas.DataFrame, palette: Palette = None, annotations: Dict[str, List[str]] = None,
			filename: Optional[Path] = None):
		if self.vertical:
			figure = plt.figure(figsize = (20, 20))
			grid = plt.GridSpec(2, 1)
			ax_timeseries: plt.Axes = figure.add_subplot(grid[0])
			ax_muller: plt.Axes = figure.add_subplot(grid[1], sharex = ax_timeseries)
		else:
			figure = plt.figure(figsize = (20, 10))
			grid = plt.GridSpec(1, 2, figure = figure)
			ax_timeseries: plt.Axes = figure.add_subplot(grid[0])
			ax_muller: plt.Axes = figure.add_subplot(grid[1], sharey = ax_timeseries)

		plotter_timeseries = TimeseriesPlot(render = self.render, style = 'nature')
		plotter_timeseries.set_scale(1)
		plotter_muller = MullerPlot(outlines = True, render = self.render, style = 'nature')
		plotter_muller.set_scale(1)

		plotter_timeseries.plot(timeseries, palette, ax = ax_timeseries)
		plotter_muller.plot(muller_df, color_palette = palette, ax = ax_muller, annotations = annotations)

		if self.vertical:
			ax_timeseries.xaxis.set_visible(False)
		else:
			ax_muller.yaxis.set_visible(False)

		plt.tight_layout()

		if filename:
			plt.savefig(filename)

	def run_from_ggmuller(self, populations: pandas.DataFrame, edges: pandas.DataFrame, palette: Optional[Dict[str, str]] = None,
			basename: Optional[Path] = None):
		""" Generates The panel object using `population` and `edges` tables."""
		from muller.dataio import GenerateMullerDataFrame
		generator = GenerateMullerDataFrame()
		muller_df = generator.run(edges, populations)

		timeseries = convert_population_to_timeseries(populations)
		return self.plot(timeseries, muller_df, palette)


def convert_population_to_timeseries(populations: pandas.DataFrame) -> pandas.DataFrame:
	""" Converts a `population` table into a `timeseries` table. """
	table = populations.pivot(index = 'Identity', columns = 'Generation', values = 'Population')
	table = table[sorted(table.columns)]

	# Need to make sure the values are within [0,1]
	table = table / 100
	return table


def _generate_image_size(size_muller: Tuple[int, int], size_lineage: Tuple[int, int], vertical: bool) -> Tuple[
	Tuple[int, int], Tuple[int, int], Tuple[int, int]]:
	muller_width, muller_height = size_muller
	lineage_width, lineage_height = size_lineage

	if vertical:
		scale = muller_width / lineage_width
		new_lineage_width = round(scale * lineage_width)
		new_lineage_height = round(scale * lineage_height)

		panel_width = new_lineage_width
		panel_height = muller_height + new_lineage_height

		position = (0, muller_height)
	else:
		scale = muller_height / lineage_height
		new_lineage_width = round(scale * lineage_width)
		new_lineage_height = round(scale * lineage_height)

		panel_width = muller_width + new_lineage_width
		panel_height = new_lineage_height
		position = (muller_width, 0)

	return (new_lineage_width, new_lineage_height), (panel_width, panel_height), position


def generate_lineage_panel(filename_muller: Path, filename_lineage: Path, filename: Path, vertical: bool = True):
	""" Combines the muller plot and lineage plots into a single figure."""

	try:
		from PIL import Image
	except ModuleNotFoundError:
		# Cannot combine the two plots.
		return None
	image_muller: Image = Image.open(filename_muller)
	image_lineage: Image = Image.open(filename_lineage)

	new_lineage_size, panel_size, position = _generate_image_size(image_muller.size, image_lineage.size, vertical)
	# Resize the lineageplot
	image_lineage = image_lineage.resize(new_lineage_size)
	# Generate the new figure
	newimage = Image.new(image_muller.mode, panel_size)
	newimage.paste(image_muller, (0, 0))
	newimage.paste(image_lineage, position)

	newimage.save(filename)


if __name__ == "__main__":
	pass
