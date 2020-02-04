import seaborn

from muller.graphics.palettes import palette_lineage


def test_apply_clade_colorscheme():
	clade = ['genotype-1', 'genotype-2', 'genotype-3']

	result = palette_lineage.apply_clade_colorscheme(clade, 'Blues')
	expected_result = dict(zip(clade, seaborn.color_palette('Blues', 3).as_hex()))
	assert result == expected_result
