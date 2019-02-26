import csv
import logging
import random
from collections import OrderedDict
from pathlib import Path
from typing import Collection, Dict, Optional, Tuple

import pandas
import seaborn

logger = logging.getLogger(__file__)


def parse_tree(edges: pandas.DataFrame) -> pandas.DataFrame:
	"""
		Determines the clade and distance from root of every leaf and node in the ggmuller edges table.
	Parameters
	----------
	edges

	Returns
	-------
	pandas.DataFrame
		- index: str
			The leaf/node genotype label
		- columns
			- 'clade': str
				The root genotype that forms the base of the clade.
			- ''distance': int
				The distance from the node/leaf to the root genotype.
	"""
	edges = edges.copy(deep = True)  # To prevent unintended alterations
	leaf_table = edges.set_index('Identity')['Parent']
	clades, iterations = zip(*[determine_clade(leaf_table, i) for i in edges['Identity'].values])
	edges['clade'] = clades
	edges['iterations'] = iterations
	edges = edges.set_index('Identity')
	return edges.sort_values(by = ['clade', 'iterations'])


def determine_clade(parents: pandas.Series, label: str) -> Tuple[str, int]:
	iteration = 0
	while True:
		iteration += 1
		index = parents[label]
		if index != 'genotype-0':
			label = index
		else:
			break
	return label, iteration


def parse_genotype_palette(paletteio: Path) -> Dict[str, str]:
	""" The file should be either a list of colors or a map of genotypes to colors."""
	palette = dict()
	with paletteio.open() as palette_file:
		reader = csv.reader(palette_file, delimiter = "\t")
		for line in reader:
			logger.debug(line)
			# Check for empty lines
			try:
				key, color, *_ = line
			except:
				continue
			if color:
				palette[key] = color
	palette['genotype-0'] = '#FFFFFF'
	return palette


def generate_random_color() -> str:
	r = random.randint(50, 200)
	g = random.randint(50, 200)
	b = random.randint(50, 200)
	color = "#{:>02X}{:>02X}{:>02X}".format(r, g, b)
	return color


def generate_clade_palette(edges_table: pandas.DataFrame) -> Dict[str, str]:
	clades = parse_tree(edges_table)
	clade_groups = clades.groupby(by = 'clade')
	genotype_colors = dict()
	color_labels = ["Greens_d", "Reds_d", "Blues_d", 'Purples_d', "Greys_d", "Oranges_d"]
	color_labels += ["winter", "autumn", "copper", "pink", "cool"]
	for base_color, (clade_root, clade) in zip(color_labels, clade_groups):
		color_palette = seaborn.color_palette(base_color, len(clade))
		clade_colors = {g: c for g, c in zip(clade.index, map(rgbtohex, color_palette))}
		genotype_colors.update(clade_colors)
	genotype_colors['genotype-0'] = '#FFFFFF'
	genotype_colors['removed'] = "#000000"
	return genotype_colors


def rgbtohex(rgb: Tuple[float, float, float]) -> str:
	if rgb[0] < 1.1:  # The values are formatted as a number between 0 and 1
		red = int(rgb[0] * 256)
		green = int(rgb[1] * 256)
		blue = int(rgb[2] * 256)
	else:
		red, green, blue = rgb
	hex_string = f"#{red:>02X}{green:>02X}{blue:>02X}"
	return hex_string


def generate_genotype_palette(genotypes: Collection, palette_filename: Optional[Path] = None) -> Dict[str, str]:
	""" Assigns a unique color to each genotype."""
	color_palette = [
		'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
		'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
		'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
		'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
	]

	if len(genotypes) >= len(color_palette):
		color_palette += [generate_random_color() for _ in genotypes]
	genotype_labels = sorted(genotypes, key = lambda s: int(s.split('-')[-1]))
	# Use an OrderedDict to help with providing the correct order for the r script.
	color_map = OrderedDict()
	color_map['genotype-0'] = "#333333"
	for label, color in zip(genotype_labels, color_palette):
		color_map[label] = color
	color_map['removed'] = '#000000'

	if palette_filename:
		custom_palette = parse_genotype_palette(palette_filename)
		color_map.update(custom_palette)

	return color_map


if __name__ == "__main__":
	pass
