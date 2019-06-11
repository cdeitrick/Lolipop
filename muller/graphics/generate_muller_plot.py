"""
	Python implementation of the Muller_plot function available from ggmuller.
"""

import math
import random
from itertools import filterfalse
from pathlib import Path
from typing import Any, Dict, List, Tuple, Union

import pandas
from matplotlib import pyplot as plt
# plt.switch_backend('agg')
from matplotlib.figure import Axes  # For autocomplete

# Make sure the svg files save the labels as actualt text.
plt.rcParams['svg.fonttype'] = 'none'
try:
	from muller.widgets import calculate_luminance
except ModuleNotFoundError:
	from ..widgets import calculate_luminance

from loguru import logger

plt.style.use('seaborn-white')

pandas.set_option('display.max_rows', 500)
pandas.set_option('display.max_columns', 500)
pandas.set_option('display.width', 250)


def unique_everseen(iterable, key = None):
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


def generate_muller_series(muller_df: pandas.DataFrame, color_palette: Dict[str, str]) -> Tuple[List[float], List[List[float]], List[str], List[str]]:
	"""
		Generates the required inputs for matplotlib to generate a mullerplot.
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

	x = list(unique_everseen(muller_df['Generation'].tolist()))
	labels = [(label if not label.endswith('a') else None) for label in genotype_order]

	groups = muller_df.groupby(by = 'Group_id')
	ys = list()
	colors = list()
	seen = set()
	for index, label in enumerate(genotype_order):
		if label in seen: continue
		seen.add(label)
		genotype_color = color_palette[label[:-1] if label.endswith('a') else label]
		try:
			next_label = genotype_order[index + 1]
		except IndexError:
			next_label = "   "  # Three spaces, so the following subscription works
		if label == next_label[:-1] and next_label[-1].isalpha():  # make sure it ends with a letter.
			# The two halves of this specific genotype are next to each other. combine them and skip the next iteration.
			seen.add(next_label)
			genotype_series_a = groups.get_group(label)['Frequency'].tolist()
			genotype_series_b = groups.get_group(next_label)['Frequency'].tolist()
			genotype_series = [i + j for i, j in zip(genotype_series_a, genotype_series_b)]
		else:
			# The genotype contains child genotypes.
			genotype_series = groups.get_group(label)['Frequency'].tolist()
		ys.append(genotype_series)
		colors.append(genotype_color)

	return x, ys, colors, labels


def get_coordinates(muller_df: pandas.DataFrame) -> Dict[str, Tuple[int, float]]:
	"""
		Attempts to find the mean point of the stacked area plot for each genotype.
	Parameters
	----------
	muller_df:pandas.DataFrame

	Returns
	-------
	points: Dict[str, Tuple[int, float]]
	A dictionary mapping each genotype to an estimate of the midpoint of the genotype's area.
	"""

	genotype_order = list()
	for index, row in muller_df.iterrows():
		genotype_label = row['Group_id']
		if genotype_label not in genotype_order:
			genotype_order.append(genotype_label)

	nonzero = muller_df[muller_df['Population'] != 0]

	points = dict()
	groups = nonzero.groupby(by = 'Group_id')

	for index, name, in enumerate(genotype_order):
		genotype_label = name[:-1] if name.endswith('a') else name
		if genotype_label in points: continue
		try:
			group = groups.get_group(name)
		except KeyError:
			continue

		min_x = min(group['Generation'].values)
		max_x = max(group['Generation'].values)
		x = (min_x + max_x) / 2
		if x not in group['Generation'].values:
			x = sorted(group['Generation'].values, key = lambda s: s - x)[0]

		gdf = nonzero[nonzero['Group_id'].isin(genotype_order[:index] + [genotype_label])]
		gdf = gdf[gdf['Generation'] == x]

		mid_y = sum(gdf['Frequency'].values)
		if mid_y < 0: mid_y = 0
		if x == 0:
			x = 0.1
		points[genotype_label] = (x, mid_y)  # + .25)

	return points


def get_font_properties(genotype_color: str) -> Dict[str, Any]:
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
	if luminance > 0.5:
		font_color = "#333333"
	else:
		font_color = "#FFFFFF"
	label_properties = {
		'size':  9,
		'color': font_color
	}
	return label_properties


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
		logger.debug(f"{x_loc}\t{y_loc}\t{closest_neighbor}\t{y_is_close}\t{y_is_almost_close}\t{x_is_close}\t{x_large}")
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


def annotate_axes(ax: Axes, points: Dict[str, Tuple[float, float]], annotations: Dict[str, List[str]], color_palette: Dict[str, str]) -> Axes:
	locations = list()
	for genotype_label, point in points.items():
		background_properties = {
			'facecolor': color_palette[genotype_label],
			'alpha':     1,
			'edgecolor': "#FFFFFF",
			'linewidth': 0.5
		}
		if genotype_label == 'genotype-0': continue
		label_properties = get_font_properties(color_palette[genotype_label])

		if locations:
			x_loc, y_loc = relocate_point(point, locations)
		else:
			x_loc, y_loc = point
		locations.append((x_loc, y_loc))
		genotype_annotations = annotations.get(genotype_label, [])
		try:
			genotype_annotations = genotype_annotations
		except IndexError:
			pass
		ax.text(
			x_loc, y_loc,
			"\n".join(genotype_annotations),
			bbox = background_properties,
			fontdict = label_properties
		)
	return ax


def generate_muller_plot(muller_df: pandas.DataFrame, color_palette: Dict[str, str],
		output_filename: Union[Path, List[Path]], annotations: Dict[str, List[str]] = None, outlines:bool = True):
	"""
		Generates a muller diagram equivilent to r's ggmuller. The primary advantage is easier annotations.
	Parameters
	----------
	muller_df: pandas.DataFrame
	color_palette: Dict[str,str]
	output_filename: Path or List[Path]
		If a list of paths is given the same graphics will be saves as each filename. This is useful if you want to save the
		plot in a number of different formats.
	annotations: Dict[str, List[str]]
		A map of genotype labels to add to the plot.
	outlines: bool; default True
		Whether to add a white outline to the series in the muller plot.
	Returns
	-------
	ax: Axes
		The axes object containing the plot.
	"""
	if not isinstance(output_filename, list):
		output_filenames = [output_filename]
	else:
		output_filenames = output_filename
	if outlines:
		edgecolor = '#FFFFFF'
	else:
		edgecolor = None

	points = get_coordinates(muller_df)

	# noinspection PyUnusedLocal,PyUnusedLocal
	fig, ax = plt.subplots(figsize = (12, 10))
	ax: Axes

	x, y, colors, labels = generate_muller_series(muller_df, color_palette)

	ax.stackplot(x, y, colors = colors, labels = labels, edgecolor = edgecolor, linewidth = 2, interpolate = True, joinstyle = 'round')
	if annotations:
		annotate_axes(ax, points, annotations, color_palette)
	else:
		plt.legend(loc = 'right', bbox_to_anchor = (1.05, 0.5, 0.1, 0), title = 'Identity')
	ax.set_facecolor("#FFFFFF")
	ax.set_xlabel("Generation", fontsize = 32)
	ax.set_ylabel("Frequency", fontsize = 32)
	for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
		label.set_fontsize(20)


	# Basic stacked area chart.
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlim(0, max(x))

	for output_filename in output_filenames:
		plt.savefig(str(output_filename))

	return ax
