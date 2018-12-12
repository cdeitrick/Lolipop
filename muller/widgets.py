import random
from typing import Dict, List, Union, Optional
import re
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


def generate_genotype_palette(genotypes: pandas.Index) -> Dict[str, str]:
	color_palette = [
		'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
		'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
		'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
		'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
	]

	if len(genotypes) >= len(color_palette):
		color_palette += [generate_random_color() for _ in genotypes]

	color_map = {i: j for i, j in zip(sorted(genotypes, key = lambda s:int(s.split('-')[-1])), color_palette)}
	color_map['genotype-0'] = "#333333"
	color_map['removed'] = '#000000'

	return color_map


def map_trajectories_to_genotype(genotypes: pandas.DataFrame) -> Dict[str, str]:
	trajectory_to_genotype = dict()
	for label, genotype in genotypes.iterrows():
		items = genotype['members'].split('|')
		for i in items:
			trajectory_to_genotype[i] = label
	return trajectory_to_genotype

