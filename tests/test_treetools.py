import pytest
import treetools
import dataio


@pytest.mark.parametrize(
	"left,right,expected",
	[
		("abcdef", "defjkl", 3),
		("abcdef", "zyx", 0),
		("jklmnop", "mno", 3)
	]
)
def test_common_values(left, right, expected):
	assert treetools.common_elements(left, right) == expected


def test_tokenize():
	expected = ["abc", "def", "jkl", "mno", "pqr"]

	result = treetools.tokenize(["abc def", "jkl", "mno pqr"])

	assert result == expected


def test_parse_tree():
	string = """
		Parent	Identity
		genotype-0	genotype-1
		genotype-1	genotype-13
		genotype-1	genotype-12
		genotype-13	genotype-9
		genotype-9	genotype-10
		genotype-0	genotype-5
		genotype-1	genotype-11
		genotype-10	genotype-7
		genotype-1	genotype-4
		genotype-13	genotype-6
		genotype-0	genotype-2
		genotype-1	genotype-3
		genotype-1	genotype-8
	"""
	table = dataio.import_table(string)

	expected_parent = {
		'genotype-1':  'genotype-1',
		'genotype-2':  'genotype-2',
		'genotype-3':  'genotype-1',
		'genotype-4':  'genotype-1',
		'genotype-5':  'genotype-5',
		'genotype-6':  'genotype-1',
		'genotype-7':  'genotype-1',
		'genotype-8':  'genotype-1',
		'genotype-9':  'genotype-1',
		'genotype-10': 'genotype-1',
		'genotype-11': 'genotype-1',
		'genotype-12': 'genotype-1',
		'genotype-13': 'genotype-1'
	}
	expected_distance = {
		'genotype-1':  1,
		'genotype-2':  1,
		'genotype-3':  2,
		'genotype-4':  2,
		'genotype-5':  1,
		'genotype-6':  3,
		'genotype-7':  5,
		'genotype-8':  2,
		'genotype-9':  3,
		'genotype-10': 4,
		'genotype-11': 2,
		'genotype-12': 2,
		'genotype-13': 2
	}

	result = treetools.parse_tree(table)

	assert result['clade'].to_dict() == expected_parent
	assert result['iterations'].to_dict() == expected_distance

def test_group_clades():
	clades = {
		'genotype-10': ['dltB '],
		'genotype-11': ['dltB Q111*'],
		'genotype-16': ['PROKKA_00173|PROKKA_00174 '],
		'genotype-6':  ['SPAR113_1988 G119C'],
		'genotype-9':  ['dltB ']
	}


	result = treetools.group_clades(clades)
	assert result == [['genotype-10', 'genotype-11', 'genotype-9'], ['genotype-16'], ['genotype-6']]