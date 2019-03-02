import pytest
import unittest.mock as mock
from dataio import file_parsers, import_table
from io import StringIO
test_labels = ["PROKKA_00139/>glyA", "gapN/>glgB", "PROKKA_00438<", "PROKKA_00512", "intergenic(62/+110)", "bglK_1<", "dnaI<","malK_1<","mepA_2<", "patA_1<"]
expected_labels = [
	"PROKKA_00139|glyA", "gapN|glgB", "PROKKA_00438", "PROKKA_00512", "intergenic", "bglK_1", "dnaI", "malK_1", "mepA_2", "patA_1"
]
@pytest.mark.parametrize(
	"test_label,expected_label",
	zip(test_labels, expected_labels)
)
def test_clean_gene_label(test_label:str, expected_label:str):

	clean_result = file_parsers._clean_gene_label(test_label)
	assert expected_label == clean_result

def test_parse_annotations():
	info_table_string = """
		Trajectory	mutation	gene
		1	A>C	PROKKA_00139/>glyA
		2	T>G	gapN/>glgB
		3	G>C	PROKKA_00438<
		4	C>T	PROKKA_00487<
		5	A>T	PROKKA_00512
		6	C>A	intergenic(62/+110)
		7	G>T	bglK_1<
		8	G>T	dnaI<
	"""
	info_table = import_table(info_table_string, index = 'Trajectory')

	genotype_members = {
		'genotype-1': "1|3|5",
		'genotype-2': "2",
		'genotype-3': "6|7|8"
	}
	expected_result = {
		'genotype-1': ["PROKKA_00139|glyA", "PROKKA_00438", "PROKKA_00512"],
		'genotype-2': ["gapN|glgB"],
		'genotype-3': ["intergenic", "bglK_1", "dnaI"]
	}

	test_result = file_parsers.parse_annotations(genotype_members, info_table)

	assert expected_result == test_result

def test_parse_genotype_palette():
	palette = """
	genotype-3	#994567
	genotype-1	#D342A1
	removed	#333311	garbage1	garbage2
	"""
	class FakePath:
		def open(self):
			p = "\n".join(i.strip() for i in palette.split('\n'))
			return StringIO(p)
	expected = {
		'genotype-1': '#D342A1',
		'genotype-3': '#994567',
		'removed': '#333311'
	}
	result = file_parsers.parse_genotype_palette(FakePath())

	assert expected == result

def test_parse_known_genotypes():
	known_genotypes = """
	1,3,4
	trajectory-7
	t5,t6
	"""
	class FakePath:
		def read_text(self):
			p = (i.strip() for i in known_genotypes.split('\n'))
			p = "\n".join(i for i in p if i)
			return p

	expected = [
		["1", "3", "4"],
		["trajectory-7"],
		["t5", "t6"]
	]
	result = file_parsers.parse_known_genotypes(FakePath())

	assert expected == result

