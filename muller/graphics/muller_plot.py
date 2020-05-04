"""
	Python implementation of the Muller_plot function available from ggmuller.
"""

import math
import random
from itertools import filterfalse
from pathlib import Path
from typing import *

import pandas
from loguru import logger
from matplotlib import pyplot as plt
# plt.switch_backend('agg')
from matplotlib.figure import Axes  # For autocomplete

# Make sure the svg files save the labels as actualt text.
plt.rcParams['svg.fonttype'] = 'none'
try:
	from muller.widgets import calculate_luminance
	from muller.graphics import Palette
except (ModuleNotFoundError, ImportError):
	from ..widgets import calculate_luminance
	from .palettes import Palette

NumericArrayType = List[Union[int,float]]
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
class MullerPlot:
	"""
		Workflow for generating a muller diagram.
		Parameters
		----------
		outlines: bool
			Whether to add an outline each genotype.
		render: bool
			Whther to generate an svg or pdf render of the diagram.
	"""

	def __init__(self, outlines: bool, render: bool, style: Optional[str] = 'default', scale: int = 1):
		self.outlines = outlines
		self.render = render
		self.dpi = 250

		self.root_genotype_name = 'genotype-0'
		self.outline_color = '#FFFFFF' if self.outlines else None  # The color of the genotype outline
		self.root_genotype_color = '#FFFFFF'  # The color of the root genotype.
		self.style = style

		# Properties for the size of the plot labels.

		# Properties for the genotype annotations that will be added to the plot
		self.genotype_annotation_outline_color = "#FFFFFF"
		self.genotype_annotation_outline_linewidth = 0.5
		self.genotype_annotation_alpha = 1

		self.genotype_annotation_font_size = 9  # Size of the text in the genotype annotation box.
		self.genotype_annotation_font_color_light = "#FFFFFF"  # Font color to use when background is dark
		self.genotype_annotation_font_color_dark = "#333333"  # Font color to use when the background is light.

		self.scale = scale
		self.reference_x = 20
		self.reference_y = 20
		self.scale_max = 3  # The maximum scale that the plots can be scalled to in order to avoid image limitations.

		# Set up the fontsizes for each labeltype
		self.label_size_axis, self.label_size_title, self.label_size_ticks = self.set_scale(scale)

		# Set the interpolation method.
		# This makes the plots look a little less jarring.
		# Available options: ‘linear’, ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’
		self.interpolation_type = 'slinear'
		# Set the number of points to generate with the interpolation method
		self.interpolation_size = 100

	def _apply_style(self, ax: plt.Axes, title: Optional[str], maximum_x: float):
		"""
			Addes labels and such to the graph.
		Parameters
		----------
		ax: THe plt.Axes object the graph was plotted on.
		maximum_x: The maximum x-value observed in the dataset.
		"""
		if self.style == 'nature':
			self._set_style_nature()
		else:
			self._set_style_default()

		ax.set_facecolor(self.root_genotype_color)

		# Add labels to the x and y axes.
		if title:
			ax.set_title(title, fontsize = self.label_size_title)
		ax.set_xlabel(self.xaxis_label, fontsize = self.label_size_axis)
		ax.set_ylabel(self.yaxis_label, fontsize = self.label_size_axis)
		# Increase the font size of the frequency and timepoint labels
		ax.tick_params(axis = 'both', labelsize = self.label_size_ticks)
		# Basic stacked area chart.

		# Remove extra ticks
		if self.style == 'nature':
			ax.set_yticks([0, 0.5, 1])
		else:
			ax.spines['right'].set_visible(False)
			ax.spines['top'].set_visible(False)

		# Make sure that the width of the diagram only shows timepoints we have data for.
		ax.set_xlim(0, maximum_x)
		ax.set_ylim(0, 1)

		return ax

	def _initialize_plot(self, ax: plt.Axes) -> plt.Axes:
		"""
			A typical plot will have about 10 to 20 genotypes and 10-20 timepoints. The current defaults work well for
			datasets of that size, but need to scale the plot size for anything larger.
		"""
		if ax is None:
			figsize_x = self.scale * 12
			figsize_y = self.scale * 10

			fig, ax = plt.subplots(figsize = (figsize_x, figsize_y))
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

	def _set_style_default(self):
		self.xaxis_label = 'Generation'
		self.yaxis_label = 'Frequency'
		self.root_genotype_color = "white"  # More readable that '#FFFFFF'

	def _set_style_nature(self):
		self.xaxis_label = 'Generation'
		self.yaxis_label = 'Individuals'
		self.root_genotype_color = 'white'

	def _set_scaling_factor(self, width_x: int, width_y: int) -> int:
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
		self.scale = max(scale_x, scale_y)
		return self.scale

	@staticmethod
	def set_scale(scale: int = 1) -> Tuple[int, int, int]:
		# Made staticmethod so that pycharm doesn't complain about object properties being defined outside of __init__()
		label_size_axis = 24 * scale
		label_size_title = 42 * scale
		label_size_ticks = 18 * scale
		return label_size_axis, label_size_title, label_size_ticks

	@staticmethod
	def get_default_palette(labels: Iterable[str]) -> Dict[str, str]:
		from . import palettes

		return palettes.generate_palette(labels)

	def save_figure(self, filename: Path):
		""" Saves the diagram in every format available in self.filetypes"""
		plt.savefig(filename, dpi = self.dpi if filename.suffix != '.svg' else None)  # Not sure if setting DPI for svgs raises an error.

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
		genotype_order = list(unique_everseen(muller_df['Group_id'].tolist()))

		x_values = list(unique_everseen(muller_df['Generation'].tolist()))
		# Extract the labels for each genotype, preserving the order they will be plotted in. If the
		# genotype was split up, only keep one of the series labels.

		# labels = [(label if not label.endswith('a') else None) for label in genotype_order]
		labels = muller_df['Identity'].unique()
		colors = [color_palette[label[:-1] if (label.endswith('a') and label not in labels) else label] for label in genotype_order]

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
		return pandas.concat(merged_table, sort = True)

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
			genotype_label = name[:-1] if (name.endswith('a') and name not in muller_df['Identity'].unique()) else name

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

	def interpolate_values(self, x:NumericArrayType, y: NumericArrayType)->Tuple[NumericArrayType, NumericArrayType]:
		""" Interpolates the values of the input arrays. """
		from scipy import interpolate

		interpolation = interpolate.interp1d(x, y, kind = self.interpolation_type)

		# Generate the interolated values
		x_value_minimum = min(x)
		x_value_maximum = max(x)
		step = (x_value_maximum - x_value_minimum) / self.interpolation_size  # The increment number for the new x-axis.
		x_interpolated = [x_value_minimum + (step * i) for i in range(self.interpolation_size)]
		y_interpolated = [interpolation(i) for i in x_interpolated]

		return x_interpolated, y_interpolated




	def plot(self, muller_df: pandas.DataFrame, filename:Path = None, color_palette: Dict[str, str] = None,
			annotations: Optional[Dict[str, List[str]]] = None,
			title: Optional[str] = None, ax: Optional[plt.Axes] = None) -> plt.Axes:
		"""
			Generates a muller diagram equivilent to r's ggmuller. The primary advantage is easier annotations.
		Parameters
		----------
		muller_df: pandas.DataFrame
		color_palette: Dict[str,str]
		filename: Path
		annotations: Dict[str, List[str]]
			A map of genotype labels to add to the plot.
		title: Optional[str]
			Applies the given title to the plot
		ax: plt.Axes
			A preexisting Axes object to draw the plot on.
		Returns
		-------
		ax: Axes
			The axes object containing the plot.
		"""

		if color_palette is None:
			color_palette = self.get_default_palette(muller_df['Identity'].unique())

		points = self.get_coordinates(muller_df)

		merged_table = self.merge_muller_table(muller_df)

		# Need to scale the graph so very large datasets are still easily distinguishable.
		# Calculate what the scale should be given the dataset shape. Any dataset smaller than the reference is automatically converted to 1.
		self.set_scale()
		# Disable the white outlines if the dataset is very large.
		if self.scale > 2:
			self.outline_color = None  # Disable the white outlines.

		ax = self._initialize_plot(ax)
		x, y, colors, labels = self.generate_muller_series(merged_table, color_palette)
		y_interpolated = list()
		for y_series in y:
			x_interpolated, yi = self.interpolate_values(x, y_series)
			y_interpolated.append(yi)

		#x_interpolated, y_interpolated = self.interpolate_values(x, y)

		ax.stackplot(
			x_interpolated, y_interpolated,
			colors = colors,
			labels = labels,
			edgecolor = self.outline_color,
			linewidth = 2,
			interpolate = True,
			joinstyle = 'round'
		)

		if annotations:
			self.add_genotype_annotations_to_plot(ax, points, annotations, color_palette)

		ax = self._apply_style(ax, title, max(x))
		if filename:
			#logger.info(f"Saving the muller plot as {filename.absolute()}")
			self.save_figure(filename)

		return ax
