import csv
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional

import pandas

logger = logging.getLogger(__file__)


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
			except ValueError:
				continue
			if color:
				palette[key] = color
	return palette


def parse_known_genotypes(known_genotypes: Path) -> List[List[str]]:
	lines = known_genotypes.read_text().split('\n')
	genotypes = [line.split(',') for line in lines]
	return genotypes


def parse_gene_map(filename: Path) -> Dict[str, str]:
	try:
		contents = filename.read_text().split('\n')
	except FileNotFoundError:
		contents = []
	alias_lines = [line.split('\t') for line in contents if line]
	alias_map = {i: j for i, j in alias_lines}

	return alias_map


def _get_column_values(table: pandas.DataFrame, column: str) -> List[str]:
	if column not in table.columns:
		if column.title() in table.columns:
			column = column.title()
		else:
			return []

	column_values = table[column].tolist()
	return column_values


def _clean_gene_label(label: str) -> str:
	gene_separator = "/>"
	if gene_separator in label:
		sublabels = [_clean_gene_label(i) for i in label.split(gene_separator)]
		clean_label = "|".join(sublabels)
	else:
		pattern = "[a-zA-Z_0-9]+"
		match = re.match(pattern, label)
		if match:
			clean_label = match.group(0)
		else:
			clean_label = ""
	return clean_label


def parse_annotations(genotype_members: pandas.Series, info: pandas.DataFrame, alias_filename: Optional[Path] = None) -> Dict[str, List[str]]:
	if alias_filename:
		alias_map = parse_gene_map(alias_filename)
	else:
		alias_map = {}

	gene_column = 'gene'
	annotation_column = 'annotation'
	amino_acid_column = "amino acid"
	genotype_annotations = dict()
	for genotype_label, members in genotype_members.items():
		trajectory_labels = members.split('|')
		trajectory_subtable = info.loc[trajectory_labels]

		gene_column_values = _get_column_values(trajectory_subtable, gene_column)
		#annotation_column_values = _get_column_values(trajectory_subtable, annotation_column)
		#amino_acid_column_values = _get_column_values(trajectory_subtable, amino_acid_column)

		gene_column_values = [_clean_gene_label(i) for i in gene_column_values]
		gene_column_values = [alias_map.get(i, i) for i in gene_column_values]
		genotype_annotations[genotype_label] = gene_column_values
	return genotype_annotations
