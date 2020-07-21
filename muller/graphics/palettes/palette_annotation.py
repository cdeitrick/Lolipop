"""
	Generates a palette using the annotations available from the input trajectories table.
"""
from typing import Dict, List

from loguru import logger


def generate_annotation_palette(genotype_annotations: Dict[str, List[str]], annotation_palette: Dict[str, str]) -> Dict[str, str]:
	"""

	Parameters
	----------
	genotype_annotations: Dict[str,List[str]]
		The annotations grouped under a given genome.
	annotation_palette: Dict[str,str]
		A map of annotations to a preset color. Genotypes with no annotations will not be assigned a color.

	Returns
	-------

	"""
	palette = dict()
	for genotype_label, genotype_annotation in genotype_annotations.items():
		candidate_colors = [annotation_palette.get(i) for i in genotype_annotation]
		candidate_colors = sorted(set(i for i in candidate_colors if i))
		if len(candidate_colors) > 1:
			logger.warning(f"Found multiple candidate colors for genotype '{genotype_label}': {candidate_colors}")
		try:
			color = candidate_colors[0]
		except IndexError:
			color = None
		palette[genotype_label] = color
	return palette
