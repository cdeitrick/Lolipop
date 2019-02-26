
import palette
import re
from io import StringIO
def test_generate_genotype_palette():
	genotypes = ['genotype-1', 'genotype-33', 'gen-5']
	expected_palette = {
		'genotype-1':  '#e6194b',
		'genotype-33': '#ffe119',
		'gen-5':       '#3cb44b',
		'genotype-0':  '#FFFFFF',
		'removed':     '#000000'
	}

	assert expected_palette == palette.generate_genotype_palette(genotypes)

def test_generate_random_color():
	color = palette.generate_random_color()
	assert re.match("#[0-9A-F]{6}", color)

def test_parse_genotype_palette():
	string = """
	genotype-7	#CC3344
	genotype-5	#123456
	genotype-11	#0057FF
	"""
	string = "\n".join([i.strip() for i in string.split('\n') if i])
	fileio = StringIO(string)
	class MockPath:
		def open(self):
			return fileio
	expected = {
		'genotype-7': "#CC3344",
		'genotype-5': "#123456",
		'genotype-11': "#0057FF"
	}
	result = palette.parse_genotype_palette(MockPath())
	assert expected == result