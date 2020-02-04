from pathlib import Path
from typing import List, Tuple

import pandas

from muller import dataio


def get_key_pairs(filename: Path) -> List[Tuple[str, str]]:
	table = pandas.read_excel(filename, sheet_name = 'genotype')
	names = table['Genotype'].values
	names = [(i, j) for j in names for i in names]
	names = [i for i in names if i[0] != i[1]]  # Don;t compare a series against itself.
	return names


def helper_test_score(scorer, method: str, left: str, right: str, filename: Path):
	""" Applies the selected check using the given filename. The `genotype` sheet and score sheet must be present."""
	table_genotype = dataio.import_table(filename, sheet_name = 'genotype', index = 'Genotype')

	if method == 'greater':
		sheetname = 'score.greater'
		actual_score = scorer.calculate_score_greater(table_genotype.loc[left], table_genotype.loc[right])
	elif method == 'fixed':
		sheetname = 'score.fixed'
		actual_score = scorer.calculate_score_above_fixed(table_genotype.loc[left], table_genotype.loc[right])
	elif method == 'derivative':
		sheetname = 'score.derivative'
		actual_score = scorer.calculate_score_derivative(table_genotype.loc[left], table_genotype.loc[right])
	elif method == 'jaccard':
		sheetname = 'score.jaccard'
		actual_score = scorer.calculate_score_area(table_genotype.loc[left], table_genotype.loc[right])
	else:
		raise ValueError
	table_scores = dataio.import_table(filename, sheet_name = sheetname, index = 'name')
	expected_score = table_scores[left][right]

	return actual_score, expected_score
