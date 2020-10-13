import re
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Union

import pandas

from muller.dataio import read_map


def read_genotype_annotations(info: pandas.DataFrame, genotype_members: Dict[str, List[str]], alias_filename: Optional[Path] = None) -> Dict[
	str, List[str]]:
	if info is not None:
		genotype_annotations = parse_genotype_annotations(
			genotype_members,
			info,
			alias_filename
		)
	else:
		genotype_annotations = {}
	return genotype_annotations



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
			clean_label = label
			#clean_label = match.group(0)
		else:
			clean_label = ""
	if clean_label.endswith('<'):
		clean_label = clean_label[:-1]
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
		# Based on a specific format.
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
	""" Maps trajectories to their corresponding mutation.
		Parameters
		----------
		info: pandas.DataFrame
			The info table. Each column label should be lowercase.
		alias_filename: Will map each gene from the infotable to an alias.
	"""
	# Make the column labels lowercase so we don't have to care about capitalization

	info.columns = [str(i).lower() for i in info.columns]
	if alias_filename:
		alias_map = read_map(alias_filename)
	else:
		alias_map = {}

	column_label_gene = 'gene'
	annotation_column = 'annotation'
	if 'annotation' not in info.columns:
		annotation_column = 'mutation'

	trajectory_annotations: Dict[str,str] = dict()
	for trajectory_label, row in info.iterrows():
		# Since the columns from the info table should already be lowercase, we don't need this method anymore.
		gene_value = row.get(column_label_gene)
		# Need to clean up the gene label.
		gene = _clean_gene_label(gene_value)
		# Rename the gene according to the alias map.
		gene = alias_map.get(gene, gene)

		annotation_value = _extract_value(row, annotation_column)
		annotation = _clean_annotation_label(annotation_value)

		result = list()
		if gene:
			result.append(gene)
		if annotation:
			result.append(annotation)
		trajectory_annotations[trajectory_label] = " ".join(result)
	return trajectory_annotations


def parse_genotype_annotations(genotype_members: Mapping[str, Union[str, List[str]]], info: pandas.DataFrame,
		alias_filename: Optional[Path] = None) -> Dict[
	str, List[str]]:
	"""
		Parses the non-numeric tables form the input table and attempts to map annotations to genotypes.
	Parameters
	----------
	genotype_members: Mapping[str, List[str]]
	info: pandas.DataFrame
		The columns from the input table that do not contain frequency measurements. The annotations will be extracted from these columns, if available:
		- `gene`
		- `annotation`

	alias_filename
	"""
	info.columns = [i.lower() for i in info.columns]
	trajectory_annotations: Dict[str,str] = extract_annotations(info, alias_filename)

	genotype_annotations = dict()
	for genotype_label, members in genotype_members.items():
		# Make sure the members are not stored as a string.
		if isinstance(members, str):
			members = members.split('|')
		# Get all annotations relevant to this genotype.
		# Want to combine the annotation for each trajectory into a single str.
		member_values = list()
		for i in members:
			trajectory_annotation = trajectory_annotations.get(i, "")
			if trajectory_annotation:
				if isinstance(trajectory_annotation, str):
					member_values.append(trajectory_annotation)
				else:
					member_values.append(" ".join(trajectory_annotation))
		# Remove missing annotations
		# Also remove excess whitespace
		member_values: List[str] = [i.strip() for i in member_values if i]
		genotype_annotations[genotype_label] = member_values
	return genotype_annotations

if __name__ == "__main__":
	pass