import re
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Union

import pandas


def read_map(contents: Union[None,str, Path]) -> Dict[str, str]:
	"""
		Parses the contents of any file used to specify a key-value pair. This is used when defining palettes and gene aliases.
	Parameters
	----------
	contents: Union[str,Path]
		Either a path to the file or the contents of the file.

	Returns
	-------
	Dict[str,str]
	"""
	if contents is None: return dict()
	lines = readfile(contents)

	filemap = dict()
	for line in lines:
		fields = line.split('\t')
		try:
			# We only care about the first two values in the line. Anything after that is ignored. May be used to add notes ot the file.
			# It's fine if the first row contains headers. They shouldn't interfere with any of the genotype/genes/etc.
			left, right, *_ = fields
		except ValueError:
			# Not enough values in line.
			continue
		left = left.strip()
		right = right.strip()
		if left and right:
			filemap[left] = right

	return filemap


def readfile(contents: Union[str, Path]) -> List[str]:
	"""
		Takes a file path or the contents of a file and attemts to split it into lines.
	Parameters
	----------
	contents: Union[str, Path]

	Returns
	-------
	List[str]
		The lines in the file.

	"""
	if not isinstance(contents, str):
		contents = contents.read_text()
	lines = contents.split('\n')  # Use this instead of splitlines sine splitlines sometimes doesn't work right.
	# Remove any whitespace to the left and right of the contents in the line. This should avaid issues with extra whitespace padding on the lines.
	clean_lines = [l.strip() for l in lines]
	# Remove empty lines
	nonempty_lines = [l for l in clean_lines if l]
	return nonempty_lines


def parse_known_genotypes(known_genotypes: Union[str, Path]) -> List[List[str]]:
	lines = readfile(known_genotypes)

	genotypes = list()
	for line in lines:
		trajectories = line.split(',')
		trajectories = [i.strip() for i in trajectories]
		genotypes.append(trajectories)
	return genotypes


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
	"""
		Removes extra stuff included in the annotation form breseq.
	"""
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


def _extract_value(row: Mapping[str, str], column: str) -> Optional[str]:
	"""
		Extracts a value from a Mapping (usually a pandas.Series object) using a few variations of the input column.
		For example, the column containing gene labels could be names 'Gene' or 'gene'. pandas.Series objects save null values as math.nan, a float.
	"""
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


def extract_annotations(info: pandas.DataFrame, alias_filename: Optional[Path] = None) -> Dict[str, str]:
	if alias_filename:
		alias_map = read_map(alias_filename)
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


def parse_genotype_annotations(genotype_members: Mapping[str, Union[str, List[str]]], info: pandas.DataFrame,
		alias_filename: Optional[Path] = None) -> Dict[
	str, List[str]]:
	trajectory_annotations = extract_annotations(info, alias_filename)

	genotype_annotations = dict()
	for genotype_label, members in genotype_members.items():
		if isinstance(members, str):
			members = members.split('|')
		member_values: List[Optional[str]] = [trajectory_annotations.get(i) for i in members]
		# Remove missing annotations
		member_values: List[str] = [i for i in member_values if i]
		genotype_annotations[genotype_label] = member_values
	return genotype_annotations
