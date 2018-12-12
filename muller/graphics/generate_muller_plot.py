"""
	Python implementation of the Muller_plot function available from ggmuller.
"""
from itertools import filterfalse
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas
from matplotlib import pyplot as plt
# plt.switch_backend('agg')
from matplotlib.figure import Axes  # For autocomplete

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

	colors = [color_palette[genotype_label[:-1] if genotype_label.endswith('a') else genotype_label] for genotype_label in genotype_order]
	labels = [(label if not label.endswith('a') else None) for label in genotype_order]
	groups = muller_df.groupby(by = 'Group_id')
	y = [groups.get_group(label)['Frequency'].tolist() for label in genotype_order]
	"""
	for genotype_label in genotype_order:
		ddf = muller_df[muller_df['Group_id'] == genotype_label]
		y.append(ddf['Frequency'].tolist())
	"""
	return x, y, colors, labels


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
		# gdf = gdf[gdf['Population'] != 50]
		gdf = gdf[gdf['Generation'] == x]

		mid_y = sum(gdf['Frequency'].values)
		if mid_y < 0: mid_y = 0
		if x == 0:
			x = 0.1
		points[genotype_label] = (x, mid_y)  # + .25)

	return points


def get_font_properties(genotype_label: str, colormap: Dict[str, str]) -> Dict[str, Any]:
	"""
		Generates font properties for each annotations. The main difference is font color,
		which is determined by the background color.
	Parameters
	----------
	genotype_label: str
	colormap: Dict[str,str]

	Returns
	-------
	fontproperties: Dict[str, Any]
	"""
	black = [
		"#FFE119", "#42D4F4", "#BFEF45", "#FABEBE", "#E6BEFF",
		"FFFAC8", "#AAFFC3", "#FFD8B1", "#FFFFFF", "#BCF60C",
		'#46f0f0'
	]
	black = [i.lower() for i in black]

	if colormap[genotype_label] in black:
		font_color = "#333333"
	else:
		font_color = "#FFFFFF"

	label_properties = {
		'size':  9,
		'color': font_color
	}
	return label_properties


def generate_muller_plot(muller_df: pandas.DataFrame, trajectory_table: Optional[pandas.DataFrame], color_palette: Dict[str, str],
		output_filename: Path,
		annotate_all: bool):
	"""
		Generates a muller diagram equivilent to r's ggmuller. The primary advantage is easier annotations.
	Parameters
	----------
	muller_df: pandas.DataFrame
	trajectory_table: pandas.DataFrame
	color_palette: Dict[str,str]
	output_filename: Path
	annotate_all:bool

	Returns
	-------
	ax: Axes
		The axes object containing the plot.
	"""
	points = get_coordinates(muller_df)

	# noinspection PyUnusedLocal,PyUnusedLocal
	fig, ax = plt.subplots(figsize = (12, 10))
	ax: Axes

	x, y, colors, labels = generate_muller_series(muller_df, color_palette)

	ax.stackplot(x, y, colors = colors, labels = labels)

	if trajectory_table is not None and 'Gene' in trajectory_table:
		groups = trajectory_table.groupby(by = 'genotype')

		for genotype_label, point in points.items():
			label_properties = get_font_properties(genotype_label, color_palette)
			x_loc, y_loc = point

			try:
				group = groups.get_group(genotype_label)
			except KeyError:
				continue
			group = group[~group['Gene'].isnull()]
			if group.empty:
				continue
			frequencies = group[[i for i in group.columns if i in muller_df['Generation']]]
			mean: pandas.Series = frequencies.mean(axis = 1)

			if not annotate_all:
				largest = mean.nlargest(3)

				# genes = [str(i) for i in group.loc[largest.index]['Gene'].tolist() if str(i) != 'nan']
				gene_df = group.loc[largest.index]['Gene']
			else:
				# genes = [str(i) for i in group['Gene'].tolist()]
				gene_df = group['Gene']
			genes = [str(i) for i in gene_df.tolist()]
			if not genes:  # No applicable genes found
				continue

			genes = [(g[:-1] if g.endswith('<') else g) for g in genes]
			plt.text(
				x_loc, y_loc,
				"\n".join(genes),
				bbox = dict(facecolor = color_palette[genotype_label], alpha = 1),
				fontdict = label_properties
			)

	ax.set_facecolor("#FFFFFF")
	ax.set_xlabel("Generation")
	ax.set_ylabel("Frequency")

	# Basic stacked area chart.

	plt.legend(loc = 'right', bbox_to_anchor = (1.05, 0.5, 0.1, 0), title = 'Identity')
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.set_xlim(0, max(x))
	ax.set_ylim(0, 1)
	plt.tight_layout()
	plt.savefig(str(output_filename))
	plt.savefig(str(output_filename.with_suffix('.pdf')))
	return ax


if __name__ == "__main__":
	pass
