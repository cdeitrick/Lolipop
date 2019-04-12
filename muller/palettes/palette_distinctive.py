from collections import OrderedDict
from typing import Dict, List

from . import colorset


def generate_distinctive_palette(genotypes: List[str]) -> Dict[str, str]:
	""" Assigns a unique color to each genotype."""
	color_palette = colorset.get_distinctive_palette(len(genotypes))

	genotype_labels = sorted(genotypes, key = lambda s: int(s.split('-')[-1]))
	# Use an OrderedDict to help with providing the correct order for the r script.
	color_map = OrderedDict()
	for label, color in zip(genotype_labels, color_palette):
		color_map[label] = color

	return color_map
