import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Union
import pandas

logger = logging.getLogger(__file__)

def openfile(contents: Union[str,Path])->str:
	if not isinstance(contents, str):
		contents = contents.read_text()
	return contents

def read_palette(paletteio: Union[str, Path]) -> Dict[str, str]:
	""" The file should be either a list of colors or a map of genotypes to colors."""
	contents = openfile(paletteio)

	palette = dict()

	for line in contents.split('\n'):
		fields = line.split('\t')
		try:
			key, color, *_ = fields
		except ValueError:
			# Not enough values in line.
			continue
		if color:
			palette[key.strip()] = color.strip()

	return palette


def parse_known_genotypes(known_genotypes: Union[str, Path]) -> List[List[str]]:
	contents = openfile(known_genotypes)

	genotypes = list()
	for line in contents.split('\n'):
		trajectories = line.split(',')
		trajectories = [i.strip() for i in trajectories]
		genotypes.append(trajectories)
	return genotypes


def parse_gene_map(filename: Path) -> Dict[str, str]:
	contents = openfile(filename)

	alias_lines = [line.split('\t') for line in contents.split('\n') if line]
	alias_map = {i.strip(): j.strip() for i, j in alias_lines}

	return alias_map


def _get_column_values(table: pandas.DataFrame, column: str) -> List[str]:
	if column not in table.columns:
		if column.title() in table.columns:
			column = column.title()
		else:
			return []

	column_values = table[column].tolist()
	column_values = [(i if isinstance(i, str) else "") for i in column_values]
	return column_values


def _clean_gene_label(label: str) -> str:
	if not isinstance(label, str):
		return ""
	gene_separator = "/>"
	if gene_separator in label:
		sublabels = [_clean_gene_label(i) for i in label.split(gene_separator)]
		clean_label = "|".join(sublabels)
	else:
		pattern = "[a-zA-Z_0-9\s]+"
		match = re.match(pattern, label)
		if match:
			clean_label = match.group(0)
		else:
			clean_label = ""
	return clean_label


def _clean_annotation_label(label: str) -> str:
	if not isinstance(label, str): return ""
	if 'intergenic' in label:
		result = 'intergenic'
	elif 'del' in label:
		result = 'del'
	else:
		result = label.split('(')[0]
		if result == "":
			result = label
	return result


def _extract_value(row: pandas.Series, column: str) -> Optional[str]:
	if column in row:
		value = row[column]
	elif column.lower() in row:
		value = row[column.lower()]
	elif column.title() in row:
		value = row[column.title()]
	else:
		value = None
	if isinstance(value, float):
		value = None
	return value


def extract_annotations(info: pandas.DataFrame, alias_filename: Optional[Path] = None) -> Dict[str, List[str]]:
	if alias_filename:
		alias_map = parse_gene_map(alias_filename)
	else:
		alias_map = {}
	gene_column = 'gene'
	annotation_column = 'annotation'

	trajectory_annotations = dict()
	for trajectory_label, row in info.iterrows():
		gene_value = _extract_value(row, gene_column)
		annotation_value = _extract_value(row, annotation_column)
		gene = _clean_gene_label(gene_value)
		gene = alias_map.get(gene, gene)
		annotation = _clean_annotation_label(annotation_value)
		value = " ".join(i for i in [gene, annotation] if i)
		trajectory_annotations[trajectory_label] = value
	return trajectory_annotations


def parse_genotype_annotations(genotype_members: pandas.Series, info: pandas.DataFrame, alias_filename: Optional[Path] = None) -> Dict[str, List[str]]:
	trajectory_annotations = extract_annotations(info, alias_filename)
	from pprint import pprint

	genotype_annotations = dict()
	for genotype_label, members in genotype_members.items():
		if isinstance(members, str):
			members = members.split('|')
		member_values = [trajectory_annotations.get(i) for i in members]
		member_values = [i for i in member_values if i]
		genotype_annotations[genotype_label] = member_values
	return genotype_annotations
