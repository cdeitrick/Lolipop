import csv
import logging
from pathlib import Path
from typing import Dict, List

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
			except:
				continue
			if color:
				palette[key] = color
	return palette


def parse_known_genotypes(known_genotypes: Path) -> List[List[str]]:
	lines = known_genotypes.read_text().split('\n')
	genotypes = [line.split(',') for line in lines]
	return genotypes
