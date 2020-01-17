from pathlib import Path
import pytest
import pandas

from muller.dataio import annotations

folder_data = Path()

@pytest.fixture
def infotable()->pandas.DataFrame:
	folder = Path.home() / "storage" / "projects" / "muller" / "manuscript" / "data" / "traverse"
	filename = folder / "traversedata.filtered.tsv"

	table = pandas.read_table(filename).set_index('Trajectory')
	infotable = table[[i for i in table if not isinstance(i, int)]]
	infotable.columns = [i.lower() for i in infotable.columns]
	return infotable

def test_parse_genotype_annotations(infotable):
	genotype_members = {
		'genotype-1':[2, 3, 12],
		'genotype-2': [4],
		'genotype-3': [9, 10, 11, 13]
	}

	expected = {
		'genotype-1': ["rpfR Y355D", "OGDH R204S", "wspD L35P"],
		'genotype-2': ['LysR-like -6:CGATGC'],
		'genotype-3': ["d49 d49", "DUF88 A209P", "MltA NC", "wspA A407V"]
	}
	result = annotations.parse_genotype_annotations(genotype_members, infotable)

	# Simple test, since it's just a list of str.
	assert result['genotype-1'] == expected['genotype-1']

	# Full test
	assert result == expected





if __name__ == "__main__":
	pass