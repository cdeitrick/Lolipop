from muller.palettes.palette_annotation import generate_annotation_palette


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
