"""
	Python implementation of the Muller_plot function available from ggmuller.
"""

import math
import random
from itertools import filterfalse
from pathlib import Path
from typing import Any, Callable, Dict, Generator, Iterable, List, Optional, Tuple

import pandas
from loguru import logger
from matplotlib import pyplot as plt
# plt.switch_backend('agg')
from matplotlib.figure import Axes  # For autocomplete

# Make sure the svg files save the labels as actualt text.
plt.rcParams['svg.fonttype'] = 'none'
try:
	from muller.widgets import calculate_luminance
except ModuleNotFoundError:
	from ..widgets import calculate_luminance

plt.style.use('seaborn-white')

pandas.set_option('display.max_rows', 500)
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 250)


###################################################################################################################################################
################################################################# Coordinates #####################################################################
###################################################################################################################################################


def unique_everseen(iterable: Iterable[Any], key: Optional[Callable[[Any], str]] = None) -> Generator[Any, None, None]:
	"""	List unique elements, preserving order. Remember all elements ever seen."
		unique_everseen('AAAABBBCCDAABBB') --> A B C D
		unique_everseen('ABBCcAD', str.lower) --> A B C D
	"""
	seen = set()
	seen_add = seen.add
	if key is None:
		for element in filterfalse(seen.__contains__, iterable):
			seen_add(element)
			yield element
	else:
		for element in iterable:
			k = key(element)
			if k not in seen:
				seen_add(k)
				yield element


def distance(left: Tuple[float, float], right: Tuple[float, float]) -> float:
	xl, yl = left
	xr, yr = right

	num = (xl - xr) ** 2 + (yl - yr) ** 2
	return math.sqrt(num)


def find_closest_point(point, points: List[Tuple[float, float]]) -> Tuple[float, float]:
	return min(points, key = lambda s: distance(s, point))


def relocate_point(point: Tuple[float, float], locations: List[Tuple[float, float]]) -> Tuple[float, float]:
	x_loc, y_loc = point
	for _ in range(10):
		closest_neighbor = find_closest_point((x_loc, y_loc), locations)
		y_is_close = math.isclose(y_loc, closest_neighbor[1], abs_tol = .1)
		y_is_almost_close = math.isclose(y_loc, closest_neighbor[1], abs_tol = .3)
		x_is_close = math.isclose(x_loc, closest_neighbor[0], abs_tol = 2)
		x_large = x_loc >= closest_neighbor[0]
		if distance(closest_neighbor, (x_loc, y_loc)) > 1: break
		if y_is_close:
			if y_loc > closest_neighbor[1]:
				y_loc += random.uniform(.1, .2)
			else:
				y_loc -= random.uniform(.05, .10)
		elif y_is_almost_close and x_is_close and x_large:
			x_loc += .5
			break
		else:
			break
	return x_loc, y_loc


###################################################################################################################################################
################################################################# Graphics #####################################################################
###################################################################################################################################################
class BaseGenerateMullerDiagram:
	"""
		Workflow for generating a muller diagram.
		Parameters
		----------
		outlines: bool
			Whether to add an outline each genotype.
		render: bool
			Whther to generate an svg or pdf render of the diagram.
	"""

	def __init__(self, outlines: bool, render: bool):
		self.outlines = outlines
		self.render = render

		self.root_genotype_name = 'genotype-0'

		self.outline_color = '#FFFFFF' if self.outlines else None  # The color of the genotype outline
		self.root_genotype_color = '#FFFFFF'  # The color of the root genotype.
		self.filetypes = ['.png']
		if self.render: self.filetypes += ['.svg']

		# Properties for the genotype annotations that will be added to the plot
		self.genotype_annotation_outline_color = "#FFFFFF"
		self.genotype_annotation_outline_linewidth = 0.5
		self.genotype_annotation_alpha = 1

		self.genotype_annotation_font_size = 9  # Size of the text in the genotype annotation box.
		self.genotype_annotation_font_color_light = "#FFFFFF"  # Font color to use when background is dark
		self.genotype_annotation_font_color_dark = "#333333"  # Font color to use when the background is light.

		self.scale_max = 3  # The maximum scale that the plots can be scalled to in order to avoid image limitations.

	def run(self, muller_df: pandas.DataFrame, color_palette: Dict[str, str], basename: Path, annotations: Optional[Dict[str, List[str]]] = None,
			title: Optional[str] = None):
		"""
			Generates a muller diagram equivilent to r's ggmuller. The primary advantage is easier annotations.
		Parameters
		----------
		muller_df: pandas.DataFrame
		color_palette: Dict[str,str]
		basename: Path
			The base filenme of the output files. The filetypes will be taken from the available suffixes.
		annotations: Dict[str, List[str]]
			A map of genotype labels to add to the plot.
		Returns
		-------
		ax: Axes
			The axes object containing the plot.
		"""
		if len(basename.suffix) == 4:
			# Probably ends with a specific extension. Remove it since the extensions are determined from self.filetypes
			basename = basename.parent / basename.stem
		points = self.get_coordinates(muller_df)

		merged_table = self.merge_muller_table(muller_df)

		# Need to scale the graph so very large datasets are still easily distinguishable.
		number_of_timepoints = len(muller_df['Generation'].unique())
		number_of_series = len(muller_df['Group_id'].unique())

		# The current values seem to work well with datasets of 10 to 20 genotypes and 15 to 20 timepoints.
		reference_x = 20
		reference_y = 20

		# Calculate what the scale should be given the dataset shape. Any dataset smaller than the reference is automatically converted to 1.
		scale_x = math.ceil(number_of_timepoints / reference_x)
		scale_y = math.ceil(number_of_series / reference_y)

		# Make sure the graph is not too big to be plotted.
		if scale_x > self.scale_max:
			logger.debug(f"Mullerplot: Reducing x scale from {scale_x} to {self.scale_max}.")
			scale_x = self.scale_max
		if scale_y > self.scale_max:
			logger.debug(f"Mullerplot: Reducing y scale from {scale_y} to {self.scale_max}.")
			scale_y = self.scale_max

		# Disable the white outlines if the dataset is very large.
		if scale_y or scale_x > 2:
			self.outline_color = None  # Disable the white outlines.

		fig, ax = self._init_plot(scale_x, scale_y)
		x, y, colors, labels = self.generate_muller_series(merged_table, color_palette)

		ax.stackplot(
			x, y,
			colors = colors,
			labels = labels,
			edgecolor = self.outline_color,
			linewidth = 2,
			interpolate = True,
			joinstyle = 'round'
		)

		if annotations:
			self.add_genotype_annotations_to_plot(ax, points, annotations, color_palette)

		ax = self._format_plot(ax, title, max(x), scale_x, scale_y)

		self.save_figure(basename)

		return ax

	@staticmethod
	def _init_plot(size_x: int, size_y: int) -> Tuple[plt.Figure, plt.Axes]:
		"""
			A typical plot will have about 10 to 20 genotypes and 10-20 timepoints. The current defaults work well for
			datasets of that size, but need to scale the plot size for anything larger.
		"""
		figsize_x = size_x * 12
		figsize_y = size_y * 10

		fig, ax = plt.subplots(figsize = (figsize_x, figsize_y))
		return fig, ax

	def _format_plot(self, ax: plt.Axes, title: Optional[str], maximum_x: float, scale_x: int, scale_y: int):
		"""
			Addes labels and such to the graph.
		Parameters
		----------
		ax: THe plt.Axes object the graph was plotted on.
		maximum_x: The maximum x-value observed in the dataset.
		scale_x, scale_y: What the plot labels should be scalled by to be easily visible on a large plot.
		"""
		scale = min(scale_x, scale_y)
		axes_label_fontsize = 32 * scale
		axes_title_fontsize = 40 * scale
		axes_label_tick_fontsize = 20 * scale

		ax.set_facecolor(self.root_genotype_color)

		# Add labels to the x and y axes.
		if not title:
			title = "Genotype Frequency"
		ax.set_title(title, fontsize = axes_title_fontsize)
		ax.set_xlabel("Generation", fontsize = axes_label_fontsize)
		ax.set_ylabel("Frequency", fontsize = axes_label_fontsize)

		# Increase the font size of the frequency and timepoint labels
		for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
			label.set_fontsize(axes_label_tick_fontsize)

		# Basic stacked area chart.
		# Remove extra ticks
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)

		# Make sure that the width of the diagram only shows timepoints we have data for.
		ax.set_xlim(0, maximum_x)

		return ax

	def save_figure(self, basename: Path):
		""" Saves the diagram in every format available in self.filetypes"""
		for suffix in self.filetypes:
			filename = str(basename) + suffix
			plt.savefig(filename, dpi = 500 if suffix != '.svg' else None)  # Not sure if setting DPI for svgs raises an error.

	# TODO: Add metadata to png images? https://matplotlib.org/3.1.1/api/backend_agg_api.html#matplotlib.backends.backend_agg.FigureCanvasAgg.print_png
	@staticmethod
	def generate_muller_series(muller_df: pandas.DataFrame, color_palette: Dict[str, str]) -> Tuple[
		List[float], List[List[float]], List[str], List[str]]:
		"""
			Generates the required inputs for matplotlib to generate a lineage.
		Parameters
		----------
		muller_df: pandas.DataFrame
		color_palette: Dict[str,str]
			Maps genotypes to a specific color.

		Returns
		-------
		x, y, colors, labels
		"""
		# TODO Add a method to highlight clades based on a gene or mutation of interest.
		genotype_order = list(unique_everseen(muller_df['Group_id'].tolist()))

		x_values = list(unique_everseen(muller_df['Generation'].tolist()))
		# Extract the labels for each genotype, preserving the order they will be plotted in. If the
		# genotype was split up, only keep one of the series labels.
		labels = [(label if not label.endswith('a') else None) for label in genotype_order]
		colors = [color_palette[label[:-1] if label.endswith('a') else label] for label in genotype_order]
		groups = muller_df.groupby(by = 'Group_id')

		# Keep track of which genotype series have already been processed.
		# Skip the second genotype series if both halves have already been merged.
		ys = list()
		for index, label in enumerate(genotype_order):
			# The genotype contains child genotypes, so it has to be kept split.
			genotype_series = groups.get_group(label)['Frequency'].tolist()
			ys.append(genotype_series)

		return x_values, ys, colors, labels

	def merge_muller_table(self, muller_df: pandas.DataFrame) -> pandas.DataFrame:
		""" The muller dataframe splits genotypes into halves so they can be plotted as a normal stacked area chart.
			We want to merge these halves if they are plotted adjacent to each other.
		"""
		# Get the order that each series appears.
		genotype_order = list(unique_everseen(muller_df['Group_id'].tolist()))
		# Each genotype should have a corresponding second-half table.
		groups = muller_df.groupby(by = 'Group_id')
		seen = set()
		merged_table = list()
		# Iterate over the population table and merge adjacent genotype subseries.
		for index, label in enumerate(genotype_order):
			if label in seen: continue
			seen.add(label)
			# Check if both of the series for this genotype are neighbors. If so, merge them.
			try:
				next_label = genotype_order[index + 1]
			except IndexError:
				next_label = "   "  # Three spaces, so the following subscription works
			# Check whether the genotype label is the same, and that the second label has an additional `a` character at the end.
			if label == next_label[:-1] and next_label[-1].isalpha():  # make sure it ends with a letter.
				# The two halves of this specific genotype are next to each other. combine them and skip the next iteration.
				seen.add(next_label)
				genotype_series = self._merge_groups(groups.get_group(label), groups.get_group(next_label))
			else:
				# The genotype contains child genotypes, so it has to be kept split.
				genotype_series = groups.get_group(label)
			merged_table.append(genotype_series)
		return pandas.concat(merged_table)

	@staticmethod
	def _merge_groups(left: pandas.DataFrame, right: pandas.DataFrame) -> pandas.DataFrame:
		"""
			Merges two genotype series which are ajacent to each other.
		"""
		left = left.set_index('Generation')
		right = right.set_index('Generation')
		frequencies = left['Frequency'].add(right['Frequency'], fill_value = 0)
		populations = left['Population'].add(right['Population'], fill_value = 0)

		df = pandas.DataFrame([frequencies, populations]).transpose()
		# Add the other columns to the dataframe
		reference = left.iloc[0]
		df['Group_id'] = reference['Group_id']
		df['Identity'] = reference['Identity']
		df['Unique_id'] = [f"{i}_{j}" for i, j in zip(df['Identity'].values, df.index)]

		df = df.reset_index()
		return df

	def add_genotype_annotations_to_plot(self, *args, **kwargs):
		""" Not neaded for a basic muller plot """
		return None

	def get_coordinates(self, *args, **kwargs):
		raise NotImplementedError


class AnnotatedMullerDiagram(BaseGenerateMullerDiagram):
	# Separates the logic used to annotate the muller diagram from the simpler logic used to actually plot it.
	def __init__(self, *args, **kwargs) -> None:
		super().__init__(*args, **kwargs)

	def add_genotype_annotations_to_plot(self, ax: Axes, points: Dict[str, Tuple[float, float]], annotations: Dict[str, List[str]],
			color_palette: Dict[str, str]) -> Axes:
		locations = list()
		for genotype_label, point in points.items():
			if genotype_label == self.root_genotype_name:
				# There's no point in adding an annotation for the root genotype.
				continue

			genotype_color: str = color_palette[genotype_label]
			genotype_annotations: List[str] = annotations.get(genotype_label, [])
			if not genotype_annotations:
				# No annotations for this genome. Don't draw anything.
				continue

			background_properties = self._get_annotation_label_background_properties(genotype_color)
			label_properties = self._get_annotation_label_font_properties(genotype_color)

			if locations:
				x_loc, y_loc = relocate_point(point, locations)
			else:
				x_loc, y_loc = point

			locations.append((x_loc, y_loc))
			ax.text(
				x_loc, y_loc,
				"\n".join(genotype_annotations),
				bbox = background_properties,
				fontdict = label_properties
			)
		return ax

	def _get_annotation_label_background_properties(self, genotype_color: str) -> Dict[str, str]:
		""" Returns a dictionary with the properties describing the genotype annotation background."""
		background_properties = {
			'facecolor': genotype_color,
			'alpha':     self.genotype_annotation_alpha,
			'edgecolor': self.genotype_annotation_outline_color,
			'linewidth': self.genotype_annotation_outline_linewidth
		}
		return background_properties

	def _get_annotation_label_font_properties(self, genotype_color: str) -> Dict[str, Any]:
		"""
			Generates font properties for each annotations. The main difference is font color,
			which is determined by the background color.
		Parameters
		----------
		genotype_color: str
			The color of the genotype this text is associated with.

		Returns
		-------
		fontproperties: Dict[str, Any]
		"""

		luminance = calculate_luminance(genotype_color)
		label_properties = {
			'size':  self.genotype_annotation_font_size,
			'color': self.genotype_annotation_font_color_dark if luminance > 0.5 else self.genotype_annotation_font_color_light
		}
		return label_properties

	def get_coordinates(self, muller_df: pandas.DataFrame) -> Dict[str, Tuple[int, float]]:
		"""
			Attempts to find the mean point of the stacked area plot for each genotype.
			This is used to determine where to draw any annotations available for the given genotype.
		Parameters
		----------
		muller_df:pandas.DataFrame

		Returns
		-------
		points: Dict[str, Tuple[int, float]]
		A dictionary mapping each genotype to an estimate of the midpoint of the genotype's area.
		"""

		# We need to keep track of when each series is plotted. There will be one or two series for each genotype, depending on if they have
		# been merged.
		genotype_order = list()
		for index, row in muller_df.iterrows():
			genotype_label = row['Group_id']
			if genotype_label not in genotype_order:
				genotype_order.append(genotype_label)

		# We don't care about when the series had a value of zero. Exclude these timepoints.
		nonzero = muller_df[muller_df['Population'] != 0]

		# Iterate over the nonzero timepoints for each series and try to find the middle of the series.
		points = dict()
		groups = nonzero.groupby(by = 'Group_id')

		for index, name, in enumerate(genotype_order):
			# muller_df splits each genotype so that it can draw them in the correct order as a stacked area chart.
			# The second series label has an additional 'a' character at the end to distinguish it from the first series for each genotype.
			genotype_label = name[:-1] if name.endswith('a') else name

			# Check if this genotype has already been assigned a location.
			if genotype_label in points: continue
			try:
				group = groups.get_group(name)
			except KeyError:
				# Just in case something weird happens.
				continue

			# Find the mean timepoint for this series.
			centroid = self.calculate_centroid(group.set_index('Generation')['Frequency'])

			# This is the centroid for the genotype series as-is.
			# This does not take into account the fact that the previously plotted genotypes will change the plotted
			# y-value of the series. So, find how much the y-value is changed due to plotting and add it to the centroid's y-value.
			previous_genotype_labels = genotype_order[:index]
			previous_genotype_table = muller_df[(muller_df['Group_id'].isin(previous_genotype_labels)) & (muller_df['Generation'] == centroid[0])]
			previous_frequencies = sum(previous_genotype_table['Frequency'].values)
			points[genotype_label] = (centroid[0], previous_frequencies + centroid[1])
		return points

	@staticmethod
	def calculate_centroid(timeseries: pandas.Series) -> Tuple[int, float]:
		""" Attempts to find the centroid for a timeseries. Used to figure out where to plot the genotype annotation"""
		mean_x = (max(timeseries.index) - min(timeseries.index)) / 2

		if mean_x not in timeseries.index:
			# Find the closest timepoint to this mean.
			mean_x = sorted(timeseries.index, key = lambda s: s - mean_x)[0]

		mean_y = timeseries.loc[mean_x]

		return mean_x, mean_y
