from pathlib import Path
from typing import *

import pandas


def read_map(contents: Union[None, str, Path]) -> Dict[str, str]:
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


