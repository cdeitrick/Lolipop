import matplotlib.pyplot as plt

plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Dict, Optional, Tuple
from muller import dataio, widgets
import math
from loguru import logger


class TimeseriesPlot:
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

		# Parameters concerning the overall plot
		self.default_color = "#333333"  # Usually only used for filtered trajectories
		self.xaxis_label = "Generation"
		self.yaxis_label = 'Frequency'
		self.markersize = 2  # The size of individual markers in the plot
		self.markertype = '-'  # The type of marker to use.
		self.linewidth = 1

		# Parameters concerning the plot legend.
		self.legend = legend
		self.legend_font_properties = {'size': 5}  # The size of the font used to label each series in the legend.
		self.legend_location = 'right'
		self.legend_title = 'Genotypes'

		self.dpi = 500
		self.scale_max = 5

	# noinspection PyTypeChecker
	def run(self, genotype_timeseries: pandas.DataFrame, genotypes: dataio.GenotypeCollection, basename: Path, palette_type: str,
			trajectory_timeseries: Optional[pandas.DataFrame]):
		"""
			Plots the clustered muller_genotypes.
		Parameters
		----------
		genotype_timeseries: pandas.DataFrame
		genotypes: dataio.Genotypes
			Has information related to how the trajectories relate to each genotype.
		basename: Optional[Path]
			Path to save the genotype plot. File extensions are determined automatically.
		palette_type: {'color_clade', 'color_unique'}
			Sets which palette to use.
		trajectory_timeseries: pandas.DataFrame
			The table with each trajectory used to generate the genotypes.
		"""
		if len(basename.suffix) == 4:
			# Probably ends with a specific extension. Remove it since the extensions are determined from self.filetypes
			basename = basename.parent / basename.stem
		trajectory_palette = genotypes.trajectory_palette(palette_type)

		# Need to scale the graph so very large datasets are still easily distinguishable.
		number_of_timepoints = len(genotype_timeseries.columns)
		number_of_series = len(genotype_timeseries)

		# The current values seem to work well with datasets of 10 to 20 genotypes and 15 to 20 timepoints.
		reference_x = 20
		reference_y = 20

		# Calculate what the scale should be given the dataset shape. Any dataset smaller than the reference is automatically converted to 1.
		scale_x = math.ceil(number_of_timepoints / reference_x)
		scale_y = math.ceil(number_of_series / reference_y)

		if scale_x > self.scale_max:
			logger.debug(f"Timeseries Plot: Reducing x scale from {scale_x} to {self.scale_max}.")
			scale_x = self.scale_max
		if scale_y > self.scale_max:
			logger.debug(f"Timeseries Plot: Reducing y scale from {scale_y} to {self.scale_max}.")
			scale_y = self.scale_max

		# Set up the plotting area. The genotype plot will cover all three columns while the trajectory plot will be relegated to the
		# Upper left. The legend will occupy the remaining space in the upper right.
		grid = plt.GridSpec(2, 3, wspace = 0.4, hspace = 0.3)

		# Plot all trajectories
		if trajectory_timeseries is not None:
			taxes: plt.Axes = plt.subplot(grid[0, 0:2])
			self.plot_timeseries(trajectory_timeseries, trajectory_palette, taxes)

			gaxes: plt.Axes = plt.subplot(grid[1, :])
		else:
			# The genotype plot can take up the full plot area.
			gaxes: plt.Axes = plt.subplot(grid[:, :])

		# Plot clustered genotypes. Should be same as above, but colored based on genotype cluster.
		# Plot the mean of each cluster

		self.plot_timeseries(genotype_timeseries, genotypes.get(palette_type), gaxes, scale = (scale_x, scale_y))
		# noinspection PyUnboundLocalVariable

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

		self.save_figure(basename)

	def _annotate_timeseries(self, axes: plt.Axes, plot_title) -> plt.Axes:
		axes.set_ylabel(self.yaxis_label)
		axes.set_title(plot_title)
		axes.set_ylim(0, 1.01)  # The maximum possible value for our data is 100% AKA 1. Leave a little room so the lines at 100% aren't obscured.

		return axes

	def plot_timeseries(self, timeseries: pandas.DataFrame, palette: Dict[str, str], ax: Optional[plt.Axes] = None, basename: Optional[Path] = None,
			scale: Tuple[int, int] = (1, 1)):
		""" Plots a generic timeseries dataframe. The plot labels are inferred from how the index labels are formatted.
			Parameters
			----------
			timeseries: pandas.DataFrame
			palette: Dict[str,str]
				Maps each series id to the proper color to use.
			ax: plt.Axes
				Specifies the plt.Axes object to use.
			basename: Optional[Path]
				The resulting figure will be saved to this filename if it is provided. The filetypes will be determined from the
				`self.filetypes` parameter.
			scale: Tuple[int,int]
				How much to adjust the size of the plot for large datasets.

		"""
		# Set up the plotting area.
		if ax is None:
			fig, ax = plt.subplots(figsize = (scale[0] * 10, scale[1] * 5))

		ax = self._annotate_timeseries(ax, 'Genotypes' if 'genotype' in timeseries.index[0] else 'Trajectories')

		numeric_columns = list(widgets.get_numeric_columns(timeseries.columns))

		for series_id, series in timeseries.iterrows():
			color = palette.get(series_id, self.default_color)
			# Make sure that the timeseries is in the right order
			# Some datasets may be out of order.
			trajectory_timeseries = sorted((column, series[column]) for column in numeric_columns)
			x_values, y_values = zip(*trajectory_timeseries)

			ax.plot(
				x_values, y_values,
				self.markertype,
				color = color,
				label = series_id,
				markersize = self.markersize,
				linewidth = self.linewidth
			)
		if basename:
			self.save_figure(basename)
		return ax

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
		for suffix in self.filetypes:
			filename = str(basename) + suffix
			plt.savefig(filename, dpi = self.dpi)
