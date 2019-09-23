from typing import Dict

import matplotlib.pyplot as plt

from muller.graphics import MullerPlot, TimeseriesPlot

plt.clf()
# plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Optional, Tuple
from muller import dataio
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
		self.legend_font_properties = {'size': 5}  # The size of the font used to label each series in the legend.
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

	def run(self, timeseries: pandas.DataFrame, genotypes: dataio.GenotypeCollection, basename: Path, palette_type: str,
			trajectory_timeseries: Optional[pandas.DataFrame]):
		"""
			Plots the clustered muller_genotypes.
		Parameters
		----------
		timeseries: pandas.DataFrame
			A dataframe where columns correspond to timepoints and rows correspond to individual series to plot.
		genotypes: dataio.Genotypes
			Has information related to how the trajectories relate to each genotype.
		basename: Optional[Path]
			Path to save the genotype plot. File extensions are determined automatically.
		palette_type: {'color_clade', 'color_unique'}
			Sets which palette to use.
		trajectory_timeseries: pandas.DataFrame
			The table with each trajectory used to generate the genotypes.
		"""

		trajectory_palette = genotypes.trajectory_palette(palette_type)

		# Need to scale the graph so very large datasets are still easily distinguishable.
		number_of_timepoints = len(timeseries.columns)
		number_of_series = len(timeseries)

		# Calculate what the scale should be given the dataset shape. Any dataset smaller than the reference is automatically converted to 1.
		scale_x, scale_y = self._get_scaling_factor(number_of_timepoints, number_of_series)
		# Create te plotting object

		# Set up the plotting area. The genotype plot will cover all three columns while the trajectory plot will be relegated to the
		# Upper left. The legend will occupy the remaining space in the upper right.
		plotter = TimeseriesPlot(scale = 1)
		grid = plt.GridSpec(2, 3, wspace = 0.4, hspace = 0.3)

		# Plot all trajectories
		taxes: plt.Axes = plt.subplot(grid[0, 0:2])
		plotter.plot(trajectory_timeseries, trajectory_palette, taxes)
		gaxes: plt.Axes = plt.subplot(grid[1, :])
		# Plot clustered genotypes. Should be same as above, but colored based on genotype cluster.
		# Plot the mean of each cluster
		plotter.plot(timeseries, genotypes.get(palette_type), gaxes)

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
		if basename:
			self.save_figure(basename)

	@staticmethod
	def get_legend_size(number_of_genotypes: int) -> Optional[Tuple[float, float]]:
		if number_of_genotypes < 20:
			bbox = (1, 2)
		elif number_of_genotypes < 30:
			bbox = (1, 1.8)
		else:
			bbox = None

		return bbox

	def save_figure(self, basename: Path):
		""" Saves the diagram in every format available in self.filetypes"""
		if len(basename.suffix) == 4:
			# Probably ends with a specific extension. Remove it since the extensions are determined from self.filetypes
			basename = basename.parent / basename.stem
		for suffix in self.filetypes:
			filename = str(basename) + suffix
			plt.savefig(filename, dpi = self.dpi)


class MullerPanel:
	""" A panel of graphics describing the lineage of a dataset."""

	def __init__(self, vertical:bool = True, render:bool = True):
		self.vertical = vertical  # The orientation of the plots.
		self.render = render

	def run(self, timeseries: pandas.DataFrame, muller_df: pandas.DataFrame, palette: Optional[Dict[str, str]] = None, basename: Optional[Path] = None):

		if self.vertical:
			figure = plt.figure(figsize = (20, 20))
			grid = plt.GridSpec(2, 1)
			ax_timeseries: plt.Axes = figure.add_subplot(grid[0])
			ax_muller: plt.Axes = figure.add_subplot(grid[1], sharex = ax_timeseries)
		else:
			figure = plt.figure(figsize = (20, 20))
			grid = plt.GridSpec(1, 2, figure = figure)
			ax_timeseries: plt.Axes = figure.add_subplot(grid[0])
			ax_muller: plt.Axes = figure.add_subplot(grid[1], sharey = ax_timeseries)

		plotter_timeseries = TimeseriesPlot(render = self.render, style = 'nature').set_scale(1)
		plotter_muller = MullerPlot(outlines = True, render = self.render, style = 'nature').set_scale(1)

		plotter_timeseries.plot(timeseries, palette, ax = ax_timeseries)
		plotter_muller.plot(muller_df, basename = basename, ax = ax_muller)

		if self.vertical:
			ax_timeseries.xaxis.set_visible(False)
		else:
			ax_muller.yaxis.set_visible(False)

		plt.tight_layout()

	def run_from_ggmuller(self, populations:pandas.DataFrame, edges:pandas.DataFrame, palette:Optional[Dict[str,str]] = None, basename: Optional[Path] = None):
		""" Generates The panel object using `population` and `edges` tables."""
		from muller.dataio import GenerateMullerDataFrame
		generator = GenerateMullerDataFrame()
		muller_df = generator.run(edges, populations)

		timeseries = convert_population_to_timeseries(populations)
		return self.run(timeseries, muller_df, palette, basename = basename)

def convert_population_to_timeseries(populations:pandas.DataFrame)->pandas.DataFrame:
	""" Converts a `population` table into a `timeseries` table. """
	table =  populations.pivot(index = 'Identity', columns = 'Generation', values = 'Population')
	table = table[sorted(table.columns)]

	# Need to make sure the values are within [0,1]
	table = table / 100
	return table

if __name__ == "__main__":
	filename_populations = "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/rosch/0.6.0-muller/M1/tables/pneumo.breseq.all.M1-matlab.ggmuller.populations.tsv"
	filename_edges = "/media/cld100/FA86364B863608A1/Users/cld100/Storage/projects/rosch/0.6.0-muller/M1/tables/pneumo.breseq.all.M1-matlab.ggmuller.edges.tsv"
	populations = pandas.read_csv(filename_populations, sep = '\t')
	edges = pandas.read_csv(filename_edges, sep = '\t')

	panel_generator = MullerPanel(vertical = True)
	panel_generator.run_from_ggmuller(populations, edges)
	plt.show()
