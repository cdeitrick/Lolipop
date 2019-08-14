import random
from typing import List, Tuple

import seaborn

distinctive_palette = [
	'#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231',
	'#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe',
	'#008080', '#e6beff', '#9a6324', '#fffac8', '#800000',
	'#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080'
]

similar_colorschemes = {
	'reds':   ["Reds", "Wistia", "Reds_r", "RdPu_r", "YlOrBr", "hot", "OrRd_r", "hot", "afmhot", "YlOrRd_r"],
	'blues':  ["Blues", "Purples", "BuPu", "PuBu", "cool"],
	'greens': ["Greens", "BuGn", "YlGn", "summer"],
	'greys':  ["Greys", "binary", "bone"],
	'browns': ["pink", "copper", "gist_heat"]
}

distinctive_colorschemes = [
	'Reds', 'Blues', 'Greens', 'Purples', 'hot', 'Greys', 'copper', 'cool'
]
sequential_colorschemes = [
	'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
	'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
	'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
sequential2_colorschemes = [
	'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
	'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
	'hot', 'afmhot', 'gist_heat', 'copper']
misc_colorschemes = [
	'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
	'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
	'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']
# For some reason seaborn doesn't like the misc colorschemes
distinctive_colorschemes = distinctive_colorschemes + sequential_colorschemes + sequential2_colorschemes# + misc_colorschemes


def load_colorscheme(string: str, size: int) -> List[str]:
	"""Attempts to load a given colorscheme."""

	try:
		palette = seaborn.color_palette(string, size)
	except ValueError:
		palette = seaborn.light_palette(string, size)

	return palette.as_hex()


def rgbtohex(rgb: Tuple[float, float, float]) -> str:
	if isinstance(rgb[0], float):  # The values are formatted as a number between 0 and 1
		red = int(rgb[0] * 256)
		green = int(rgb[1] * 256)
		blue = int(rgb[2] * 256)
	else:
		red, green, blue = rgb
	hex_string = f"#{red:>02X}{green:>02X}{blue:>02X}"
	return hex_string


def random_color(lower: int = 50, upper: int = 230) -> str:
	r = random.randint(lower, upper)
	g = random.randint(lower, upper)
	b = random.randint(lower, upper)
	color = "#{:>02X}{:>02X}{:>02X}".format(r, g, b)
	return color


def get_distinctive_palette(number: int = None) -> List[str]:
	""" A set of selected colors meant to be as clearly distinct from each other as possible."""
	if number:
		return distinctive_palette + [random_color() for _ in range(number)]
	return distinctive_palette
