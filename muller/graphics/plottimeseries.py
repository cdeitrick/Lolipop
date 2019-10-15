import matplotlib.pyplot as plt

plt.clf()
# plt.switch_backend('agg')
import pandas
from pathlib import Path
from typing import Dict, Optional
from muller import widgets
from loguru import logger
from muller.graphics.palettes import palette_distinctive


class TimeseriesPlot:
	def __init__(self, render: bool = True, legend: bool = True, scale: int = 1, style: Optional[str] = None):
		# Set the default aspect ratio
		self.length_x = 12
		self.length_y = 10
		self.scale = scale

		self.dpi = 250
		self.filetypes = ['.png']
		if render:
			self.filetypes += ['.svg']

		# Parameters concerning the overall plot
		self.default_color = "#333333"
		self.style = style
		if self.style == 'nature':
			self._set_style_nature()
		else:
			self._set_style_default()

		# Parameters concerning the plot legend.
		self.legend = legend
		self.legend_font_properties = {'size': 24}  # The size of the font used to label each series in the legend.
		self.legend_location = 'right'
		self.legend_title = 'Genotypes'

	def set_scale(self, scale:int = 1)->'TimeseriesPlot':
		self.scale = scale
		self.label_size_axis = 24 * self.scale
		self.label_size_title = 42 * self.scale
		self.label_size_ticks = 18 * self.scale
		return self
	def _set_style_default(self):
		# Parameters concerning the overall plot
		self.xaxis_label = "Generation"
		self.yaxis_label = 'Frequency'
		self.background_color = 'white'

		self.markersize = 2  # The size of individual markers in the plot
		self.markertype = 'o'  # The type of marker to use.

		self.linestyle = 'solid'
		self.linewidth = 12

	def _set_style_nature(self):
		""" Configures `TimeseriesPlot` to use a style similar to that in the yeast nature paper. """
		self.xaxis_label = "Generation"
		self.yaxis_label = "Frequency"
		self.background_color = "white"

		self.markertype = 'o'
		self.markersize = 12

		self.linestyle = 'solid'
		self.linewidth = 3

	def _apply_style(self, axes: plt.Axes, plot_title, xmax: Optional[int] = None) -> plt.Axes:
		# Parameters concerning subfeatures of the plot.
		axes.set_ylabel(self.yaxis_label, fontsize = self.label_size_axis)
		axes.set_xlabel(self.xaxis_label, fontsize = self.label_size_axis)
		axes.set_facecolor(self.background_color)
		axes.set_title(plot_title, fontsize = self.label_size_title)
		axes.set_ylim(0, 1.01)  # The maximum possible value for our data is 100% AKA 1. Leave a little room so the lines at 100% aren't obscured.
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

	def plot(self, timeseries: pandas.DataFrame, palette: Optional[Dict[str, str]] = None, ax: Optional[plt.Axes] = None,
			basename: Optional[Path] = None) -> plt.Axes:
		""" Plots a generic timeseries dataframe. The plot labels are inferred from how the index labels are formatted.
			Parameters
			----------
			timeseries: pandas.DataFrame
				A dataframe where each column is a timepoint and each row is a specific series to plot.
			palette: Dict[str,str]
				Maps each series id to the proper color to use.
			ax: Optional[plt.Axes]
				Specifies the plt.Axes object to use.
			basename: Optional[Path]
				The resulting figure will be saved to this filename if it is provided. The filetypes will be determined from the
				`self.filetypes` parameter.
		"""
		# Set up the plotting area.
		self.set_scale()
		if palette is None:
			palette = self.get_palette(timeseries)

		ax = self._initialize_plot(ax)
		plot_title = 'Genotypes' if 'genotype' in timeseries.index[0] else 'Trajectories'

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
				marker = self.markertype,
				markersize = self.markersize,
				linewidth = self.linewidth,
				linestyle = self.linestyle
			)
		ax = self._apply_style(ax, plot_title, max(timeseries.columns))
		ax.set_xlim(0, max(timeseries.columns))
		if self.legend and False:
			legend = ax.legend(
				loc = self.legend_location,
				prop = self.legend_font_properties,
				title = self.legend_title
			)
			legend.get_title().set_fontsize(str(self.legend_font_properties['size']))

		if basename:
			self.save_figure(basename)
		return ax

	def save_figure(self, basename: Path):
		""" Saves the diagram in every format available in self.filetypes"""
		if len(basename.suffix) == 4:
			# Probably ends with a specific extension. Remove it since the extensions are determined from self.filetypes
			basename = basename.parent / basename.stem
		for suffix in self.filetypes:
			filename = str(basename) + suffix
			plt.savefig(filename, dpi = self.dpi)


if __name__ == "__main__":


	filename = Path("/home/cld100/Documents/github/muller_diagrams/tests/data/tables/real.nature12344-s2.BYB1-G07.xlsx")
	datatable = pandas.read_excel(filename, sheet_name = 'genotype').set_index('Genotype')
	datatable.columns = [int(i) for i in datatable.columns]

	data = [
		pandas.Series([0.0, 0.1, 0.15, 0.2, 0.3, 0.5, 0.8, 0.9, 1.0, 1.0], name = 'selected'), # fixed
		pandas.Series([0.0, 0.05, 0.15, 0.1, 0.15, 0.2, 0.1, 0.1, 0.15,0.1], name = 'noise'), # low-freq
		pandas.Series([0.0, 0.2, 0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], name = 'competing'), # competing
		pandas.Series([0.0, 0.0, 0.1, 0.15, 0.25, 0.45, 0.75, 0.95, 0.95, 1.0], name = 'linked')
	]
	datatable = pandas.DataFrame(data)
	print(datatable.to_string())

	naturepalette = {
		'genotype-aqua':   'black',
		'genotype-red':    '#2ebebf',
		'genotype-orange': '#a74b9c',
		'genotype-green':  '#497fc1',
		'genotype-sienna': '#28864d',
		'genotype-gold':   '#c0c23e'
	}

	palette = {
		'selected': 'black',
		'noise': 'red',
		'competing': 'green',
		'linked': 'sienna'
	}
	plotter = TimeseriesPlot(legend = True)
	plotter.plot(datatable, palette = palette)
	plt.tight_layout()
	plt.savefig('example.filtering.png')
