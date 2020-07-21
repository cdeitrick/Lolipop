import re
from collections import OrderedDict
from typing import Dict, List

import matplotlib.colors

from . import colorset

available_colors = matplotlib.colors.get_named_colors_mapping()

def generate_distinctive_palette(genotypes: List[str]) -> Dict[str, str]:
	""" Assigns a unique color to each genotype."""
	# TODO: If the genotype label specifies a color, use that instead.
	color_palette = colorset.get_distinctive_palette(len(genotypes))
	if re.search("genotype-[\d]+", genotypes[0]):
		genotype_labels = sorted(genotypes, key = lambda s: int(s.split('-')[-1]))
	else:
		genotype_labels = genotypes
	# Use an OrderedDict to help with providing the correct order for the r script.
	color_map = OrderedDict()
	for label, color in zip(genotype_labels, color_palette):
		# If the genotype name specifies a color, use that.
		genotype_name = label.split('-')[1]
		if genotype_name in available_colors:
			color_map[label] = available_colors[genotype_name]
		else:
			color_map[label] = color

	return color_map
