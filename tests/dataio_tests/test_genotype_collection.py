import pytest

from muller import dataio


@pytest.fixture
def small_genotypes() -> dataio.GenotypeCollection:
	a = {
		'label':        'genotype-1',
		'color_clade':  '#F00000',
		'color_unique': '#FF0000',
		'color_custom': None,
		'annotations':  ['D64D(GAC>GAT)', 'P39P(CCC>CCT)'],
		'members':      ['14', '15']
	}
	b = {
		'label':        'genotype-2',
		'color_clade':  '#00F000',
		'color_unique': '#00FF00',
		'color_custom': '#333333',
		'annotations':  ["S48R(AGC>CGC)"],
		'members':      ['abc']
	}
	c = {
		'label':        'genotype-3',
		'color_clade':  '#0000F0',
		'color_unique': '#0000FF',
		'color_custom': None,
		'annotations':  ["intergenic(+12/-169)", "G317G(GGT>GGG)", "C270G(TGC>GGC)"],
		'members':      ['1', '13', '5']
	}
	gc = dataio.GenotypeCollection()
	gc['genotype-1'] = dataio.Genotype(**a)
	gc['genotype-2'] = dataio.Genotype(**b)
	gc['genotype-3'] = dataio.Genotype(**c)
	return gc


def test_palette_generator_clade(small_genotypes):
	expected = {
		'genotype-1': '#F00000',
		'genotype-2': '#333333',
		'genotype-3': '#0000F0'
	}
	assert small_genotypes.get('color_clade') == expected


def test_palette_generator_unique(small_genotypes):
	expected = {
		'genotype-1': '#FF0000',
		'genotype-2': '#333333',
		'genotype-3': '#0000FF'
	}
	assert small_genotypes.get('color_unique') == expected


def test_palette_generator_trajectory(small_genotypes):
	expected = {
		'14':  '#FF0000',
		'15':  '#FF0000',
		'abc': '#333333',
		'1':   '#0000FF',
		'13':  '#0000FF',
		'5':   '#0000FF'
	}
	assert small_genotypes.trajectory_palette('color_unique') == expected


def test_get_method_works(small_genotypes):
	expected = {
		'genotype-1': ['14', '15'],
		'genotype-2': ['abc'],
		'genotype-3': ['1', '13', '5']
	}
	assert small_genotypes.get('members') == expected


def test_get_annotation_map(small_genotypes):
	expected = {
		'D64D(GAC>GAT)':        ['genotype-1'],
		'P39P(CCC>CCT)':        ['genotype-1'],
		"S48R(AGC>CGC)":        ['genotype-2'],
		"intergenic(+12/-169)": ['genotype-3'],
		"G317G(GGT>GGG)":       ['genotype-3'],
		"C270G(TGC>GGC)":       ['genotype-3']
	}

	assert small_genotypes.get_annotation_map() == expected


def test_get_genotype_from_annotation(small_genotypes):
	assert small_genotypes.get_genotype_from_annotation('P39P(CCC>CCT)') == ['genotype-1']
	assert small_genotypes.get_genotype_from_annotation("missing") == []
