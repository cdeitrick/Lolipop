import matplotlib.pyplot as plt

plt.clf()
# plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Dict, Optional, List, Tuple, Union
from muller import widgets
from muller.graphics.palettes import palette_distinctive, Palette
from loguru import logger


class TimeseriesPlot:
	def __init__(self, render: bool = True, legend: bool = True, scale: int = 1, style: Optional[str] = None):
		# Set the default aspect ratio
		self.length_x = 12
		self.length_y = 10
		self.scale = scale
		if self.scale < 1:
			self.scale = 1

		self.dpi = 250

		# Parameters concerning the overall plot
		self.default_color = "#333333"
		self.style = style
		if self.style == 'nature':
			self._set_style_nature()
		else:
			self._set_style_default()

		# Parameters concerning the plot legend.
		self.legend = legend
		self.legend_font_properties = {
			'size': 12 * self.scale
		}  # The size of the font used to label each series in the legend.
		self.legend_location = 'right'
		self.legend_title = 'Genotypes'

		# Set up the fontsizes for each labeltype
		self.label_size_axis, self.label_size_title, self.label_size_ticks = self.set_scale(scale)

	@staticmethod
	def set_scale(scale: int = 1) -> Tuple[int, int, int]:
		# Made staticmethod so that pycharm doesn't complain about object properties being defined outside of __init__()
		label_size_axis = 24 * scale
		label_size_title = 42 * scale
		label_size_ticks = 18 * scale
		return label_size_axis, label_size_title, label_size_ticks

	def _set_style_default(self) -> None:
		# Parameters concerning the overall plot
		self.xaxis_label = "Generation"
		self.yaxis_label = 'Frequency'
		self.background_color = 'white'

		self.markersize = 4 * self.scale  # The size of individual markers in the plot
		self.markertype = 'o'  # The type of marker to use.

		self.linestyle = 'solid'
		self.linewidth = 2 * self.scale

	def _set_style_nature(self):
		""" Configures `TimeseriesPlot` to use a style similar to that in the yeast nature paper. """
		self.xaxis_label = "Generation"
		self.yaxis_label = "Frequency"
		self.background_color = "white"

		self.markertype = 'o'
		self.markersize = 12 * self.scale

		self.linestyle = 'solid'
		self.linewidth = 3 * self.scale

	def _apply_style(self, axes: plt.Axes, plot_title, xmax: Optional[int] = None) -> plt.Axes:
		# Parameters concerning subfeatures of the plot.
		axes.set_ylabel(self.yaxis_label, fontsize = self.label_size_axis)
		axes.set_xlabel(self.xaxis_label, fontsize = self.label_size_axis)
		axes.set_facecolor(self.background_color)
		axes.set_title(plot_title, fontsize = self.label_size_title)
		axes.set_ylim(0,
			1.01)  # The maximum possible value for our data is 100% AKA 1. Leave a little room so the lines at 100% aren't obscured.
		if xmax:
			axes.set_xlim(0, xmax + 1)  # SMall offset so the final point isn't sut off

		if self.style == 'nature':
			axes.set_yticks([0, 0.5, 1])

		axes.tick_params(axis = 'both', labelsize = self.label_size_ticks)

		return axes

	def _initialize_plot(self, ax: Optional[plt.axes]) -> plt.Axes:
		if ax is None:
			fig, ax = plt.subplots(figsize = (self.length_x * self.scale, self.length_y * self.scale))
		return ax

	@staticmethod
	def get_palette(table: pandas.DataFrame) -> Dict[str, str]:
		return palette_distinctive.generate_distinctive_palette(table.index)

	def plot_multiple(self, timeseries: pandas.DataFrame, palettes: List[Palette], ax: Optional[plt.Axes] = None,
			filenames: Optional[List[Path]] = None):
		for palette in palettes:
			fnames = filenames[palette.name]
			self.plot(timeseries, palette, ax, fnames)

	def plot(self, timeseries: pandas.DataFrame, palette: Union[Dict[str, str], Palette] = None,
			ax: Optional[plt.Axes] = None,
			filename: Optional[Path] = None) -> plt.Axes:
		""" Plots a generic timeseries dataframe. The plot labels are inferred from how the index labels are formatted.
			Parameters
			----------
			timeseries: pandas.DataFrame
				A dataframe where each column is a timepoint and each row is a specific series to plot.
			palette: Dict[str,str]
				Maps each series id to the proper color to use.
			ax: Optional[plt.Axes]
				Specifies the plt.Axes object to use.
			filename: Optional[Path]
				The resulting figure will be saved to this filename if it is provided.
		"""
		# Set up the plotting area.
		if palette is None: palette = {}
		self.set_scale()

		ax = self._initialize_plot(ax)
		try:
			plot_title = 'Genotypes' if 'genotype' in timeseries.index[0] else 'Trajectories'
		except TypeError:
			message = f"Could not iterate over the first element of the index. Is the table indexed by genotype?"
			logger.debug(message)
			plot_title = "Timeseries"
		numeric_columns = list(widgets.get_numeric_columns(timeseries.columns))
		timeseries = timeseries[numeric_columns]

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
				marker = self.markertype,
				markersize = self.markersize,
				linewidth = self.linewidth,
				linestyle = self.linestyle
			)
		ax = self._apply_style(ax, plot_title, max(int(i) for i in timeseries.columns))
		ax.set_xlim(0, max(timeseries.columns))
		if self.legend and False:
			legend = ax.legend(
				loc = self.legend_location,
				prop = self.legend_font_properties,
				title = self.legend_title
			)
			legend.get_title().set_fontsize(str(self.legend_font_properties['size']))

		if filename:
			self.save_figure(filename)
		return ax

	def save_figure(self, filename: Path):
		""" Saves the diagram in every format available in self.filetypes"""
		plt.savefig(filename, dpi = self.dpi)


if __name__ == "__main__":
	pass
