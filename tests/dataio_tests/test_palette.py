import re

import pandas
import pytest

from muller import treetools
from muller.graphics import palettes
from muller.dataio import import_table
from muller.graphics.palettes.palette_annotation import generate_annotation_palette
from muller.graphics.palettes import palette_lineage
import seaborn


@pytest.fixture
def edges_table() -> pandas.Series:
	data = """
		Parent	Identity
		genotype-0	genotype-2
		genotype-0	genotype-1
		genotype-1	genotype-5
		genotype-2	genotype-6
		genotype-1	genotype-4
		genotype-4	genotype-12
		genotype-6	genotype-3
		genotype-2	genotype-9
		genotype-4	genotype-14
		genotype-0	genotype-10
		genotype-10	genotype-15
		genotype-10	genotype-16
		genotype-0	genotype-11
		genotype-2	genotype-8
		genotype-2	genotype-7
		genotype-6	genotype-13
		genotype-4	genotype-17
		genotype-4	genotype-19
		genotype-14	genotype-18
	"""
	table = import_table(data, index = 'Identity')
	table = table['Parent']
	return table


@pytest.mark.parametrize("label,expected",
	[
		("genotype-8", "genotype-2"),
		("genotype-18", "genotype-1"),
		("genotype-2", "genotype-2"),
		("genotype-5", "genotype-1"),
		("genotype-9", "genotype-2"),
		("genotype-19", "genotype-1"),
		("genotype-13", "genotype-2")
	]
)
def test_determine_clade(edges_table, label, expected):
	result, iterations = treetools.determine_clade(edges_table, label)
	assert expected == result


def test_generate_random_color():
	color = palettes.random_color()
	assert re.match("#[0-9A-F]{6}", color)


@pytest.mark.parametrize("rgb,expected",
	[
		((247, 252, 253), "#f7fcfd"),
		((229, 245, 249), "#e5f5f9"),
		((204, 236, 230), "#ccece6"),
		((153, 216, 201), "#99d8c9"),
		((102, 194, 164), "#66c2a4"),
		((65, 174, 118), "#41ae76"),
		((35, 139, 69), "#238b45"),
		((0, 109, 44), "#006d2c"),
		((0, 68, 27), "#00441b")
	]
)
def test_rgbtohex(rgb, expected):
	result = palettes.rgbtohex(rgb)
	assert expected.lower() == result.lower()


def test_generate_annotation_palette():
	annotations = {
		'genotype-1': ['gene1', 'gene7'],
		'genotype-7': ['gene5'],
		'genotype-4': []
	}
	palette = {
		'gene1': '#444444',
		'gene7': '#555555',
		'gene5': '#666666'
	}

	expected_palette = {
		'genotype-1': '#444444',
		'genotype-7': '#666666',
		'genotype-4': None
	}
	result = generate_annotation_palette(annotations, palette)

	assert result == expected_palette


def test_apply_clade_colorscheme():
	clade = ['genotype-1', 'genotype-2', 'genotype-3']

	result = palette_lineage.apply_clade_colorscheme(clade, 'Blues')
	expected_result = dict(zip(clade, seaborn.color_palette('Blues', 3).as_hex()))
	assert result == expected_result
