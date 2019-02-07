import random
import re
from collections import OrderedDict
from typing import Collection, Dict, List

import pandas

NUMERIC_REGEX = re.compile("^.?(?P<number>[\d]+)")


def generate_random_color() -> str:
	r = random.randint(50, 200)
	g = random.randint(50, 200)
	b = random.randint(50, 200)
	color = "#{:>02X}{:>02X}{:>02X}".format(r, g, b)
	return color


def get_numeric_columns(columns: List[str]) -> List[str]:
	numeric_columns = list()
	for column in columns:
		if isinstance(column, str):
			match = NUMERIC_REGEX.search(column)
			if match:
				col = match.groupdict()['number']
			else:
				continue
		else:
			col = column

		try:
			int(col)
		except (ValueError, TypeError):
			continue
		numeric_columns.append(column)
	return numeric_columns


def generate_genotype_palette(genotypes: Collection) -> Dict[str, str]:
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

	return color_map


def map_trajectories_to_genotype(genotype_members: pandas.Series) -> Dict[str, str]:
	trajectory_to_genotype = dict()
	for genotype_label, members in genotype_members.items():
		for member in members.split('|'):
			trajectory_to_genotype[member] = genotype_label
	return trajectory_to_genotype
