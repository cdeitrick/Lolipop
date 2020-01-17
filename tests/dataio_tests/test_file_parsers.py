import pandas
import pytest

from muller import dataio
from muller.dataio import import_file, import_table, annotations


@pytest.fixture
def genotypeannotations() -> pandas.DataFrame:
	string = """
	Trajectory	Chromosome	Position	Class	Mutation	Gene	Annotation	Class	Amino	Description
	1	1	36,414	SNP	G>A	speA	C127Y(TGT>TAT)	Non	C->Y	Arginine decarboxylase
	2	1	138,043	SNP	C>T	rlmCD_1	H441H(CAC>CAT)	Syn	H->H	23S rRNA (uracilC(5))methyltransferase RlmCD
	3	1	165,470	SNP	C>A	PROKKA_00173/>PROKKA_00174	intergenic(+174/91)	NC	na	Relaxase/Mobilisation nuclease domain protein/hypothetical protein
	4	1	234,888	SNP	C>A	gapN	A95E(GCA>GAA)	Non	A->E	NADPdependent glyceraldehyde3phosphate dehydrogenase
	"""
	return dataio.import_table(string, index = 'Trajectory')


def test_extract_annotations(genotypeannotations):
	expected = {
		'1': 'speA C127Y',
		'2': 'rlmCD_1 H441H',
		'3': "PROKKA_00173|PROKKA_00174 intergenic",
		'4': "gapN A95E"
	}

	result = annotations.extract_annotations(genotypeannotations)

	assert result == expected


@pytest.mark.parametrize(
	"test_label,expected_label",
	[
		("PROKKA_00139/>glyA", "PROKKA_00139|glyA"),
		("gapN/>glgB", "gapN|glgB"),
		("PROKKA_00438<", "PROKKA_00438"),
		("intergenic(62/+110)", "intergenic(62/+110)"),
		("bglK_1<", "bglK_1")
	]
)
def test_clean_gene_label(test_label: str, expected_label: str):
	clean_result = annotations._clean_gene_label(test_label)
	assert clean_result == expected_label


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
		'genotype-3': ["intergenic(62/+110)", "bglK_1", "dnaI"]
	}

	test_result = annotations.parse_genotype_annotations(genotype_members, info_table)

	assert expected_result == test_result


def test_read_map():
	expected = {
		'genotype-1': '#D342A1',
		'genotype-3': '#994567',
		'removed':    '#333311'
	}
	palette = """
	genotype-3	#994567
	genotype-1	#D342A1
	removed	#333311	garbage1	garbage2
	"""
	result = import_file.read_map(palette)

	assert expected == result

	# Test a file with garbage
	palette_with_garbage = """
	genotype-3	#994567	garbage
	shortline
		genotype-1	#D342A1	extra whitespace
	removed	#333311
	"""
	result = import_file.read_map(palette_with_garbage)
	assert result == expected


def test_parse_known_genotypes():
	known_genotypes = """
	1,3,4
	trajectory-7
	t5,t6
	"""

	expected = [
		["1", "3", "4"],
		["trajectory-7"],
		["t5", "t6"]
	]
	result = import_file.parse_known_genotypes(known_genotypes)

	assert expected == result
