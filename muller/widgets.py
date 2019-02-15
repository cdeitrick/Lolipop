import csv
import random
import re
from collections import OrderedDict
from pathlib import Path
from typing import Collection, Dict, List, Optional
import logging
logger = logging.getLogger(__file__)
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


def parse_genotype_palette(paletteio: Path) -> Dict[str, str]:
	""" The file should be either a list of colors or a map of genotypes to colors."""
	palette = dict()
	with paletteio.open() as palette_file:
		reader = csv.reader(palette_file, delimiter = "\t")
		for line in reader:
			logger.debug(line)
			#Check for empty lines
			try:
				key, color = line
			except:
				continue
			palette[key] = color

	return palette


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


def map_trajectories_to_genotype(genotype_members: pandas.Series) -> Dict[str, str]:
	trajectory_to_genotype = dict()
	for genotype_label, members in genotype_members.items():
		for member in members.split('|'):
			trajectory_to_genotype[member] = genotype_label
	return trajectory_to_genotype


def get_valid_points(left: pandas.Series, right: pandas.Series, dlimit: float, flimit: Optional[float] = None,
		inner: bool = False) -> pandas.DataFrame:
	"""
		Filters out timepoints that do not satisfy the detection criteria.
	Parameters
	----------
	left: pandas.Series
	right: pandas.Series
	dlimit: float
		Removes the timepoint if both points do no exceed this value.
	flimit: float
		If given, all points exceeding this value are treated as "undetected".
	inner: bool; default False
		If `true`, both points must exceed the detection limit for a given timepoint to be considered valid.

	Returns
	-------

	"""
	df = pandas.concat([left, right], axis = 1)
	if flimit is not None:
		# Mask the values so they are considered below the detection limit.
		# Assign a value of -1 so that they will be excluded wven if the detection limit is 0.
		left = left.mask(lambda s: s > flimit, -1)
		right = right.mask(lambda s: s > flimit, -1)
	df.columns = ['left', 'right']  # Prevent an IndexError when `left` and `right` refer to the same series.
	if inner:
		at_least_one_detected = (left > dlimit) & (right > dlimit)
	else:
		at_least_one_detected = (left > dlimit) | (right > dlimit)
	# Remove indicies where the series value falls below the detection limit. This should include the masked fixed values.
	at_least_one_detected_reduced = at_least_one_detected[at_least_one_detected]
	if at_least_one_detected_reduced.empty:
		# There are no shared timepoints between the series. Assign index_min and index_max to the same number, which will result in an empty dataframe.
		position_index_min = position_index_max = 0
	else:
		try:
			position_index_min_value = min(at_least_one_detected_reduced.index)
			position_index_max_value = max(at_least_one_detected_reduced.index)
		except TypeError:
			# The indicies are str and we can't use min() or max(). Assume the indicies are already sorted.
			position_index_min_value = at_least_one_detected_reduced.index[0]
			position_index_max_value = at_least_one_detected_reduced.index[-1]
		position_index_min = df.index.get_loc(position_index_min_value)
		position_index_max = df.index.get_loc(position_index_max_value)
		# Since we want to include the last index, increment position_index_max by one.
		position_index_max += 1
	result = df.iloc[position_index_min:position_index_max]
	return result


get_detected_points = get_valid_points


def format_linkage_matrix(Z) -> pandas.DataFrame:
	Z = pandas.DataFrame(Z, columns = ["left", "right", "distance", "observations"])
	# Z.index = pandas.Index([i + len(squaremap.index) for i in Z.index], name = "clusterLabel")
	Z['left'] = Z['left'].astype(int)
	Z['right'] = Z['right'].astype(int)
	Z['observations'] = Z['observations'].astype(int)
	return Z


def format_inconsistency_matrix(R) -> pandas.DataFrame:
	inconsistency_table = pandas.DataFrame(R, columns = ['mean', 'std', 'observations', 'statistic'])
	inconsistency_table['observations'] = inconsistency_table['observations'].astype(int)
	return inconsistency_table
