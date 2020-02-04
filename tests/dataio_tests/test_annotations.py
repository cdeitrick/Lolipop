from pathlib import Path
import pytest
import pandas

from muller.dataio import annotations

folder_data = Path()

@pytest.fixture
def infotable()->pandas.DataFrame:
	folder = Path.home() / "storage" / "projects" / "muller" / "datasets" / "traverse"
	filename = folder / "traversedata.filtered.tsv"

	table = pandas.read_table(filename).set_index('Trajectory')
	infotable = table[[i for i in table if not isinstance(i, int)]]
	infotable.columns = [i.lower() for i in infotable.columns]
	return infotable

def test_parse_genotype_annotations(infotable):
	genotype_members = {
		'genotype-1':["M2", "M3", "M12"],
		'genotype-2': ["M4"],
		'genotype-3': ["M9", "M10", "M11", "M13"]
	}

	expected = {
		'genotype-1': ["monoxygenase E481D", "OGDH R204S", "wspD L35P"],
		'genotype-2': ['LysR-like -6:CGATGC'],
		'genotype-3': ["succinate dehydrogenase G147G", "DUF88 A209P", "MltA -10bp", "wspA A407V"]
	}
	result = annotations.parse_genotype_annotations(genotype_members, infotable)

	# Simple test, since it's just a list of str.
	assert result['genotype-1'] == expected['genotype-1']
	assert result['genotype-2'] == expected['genotype-2']
	assert result['genotype-3'] == expected['genotype-3']
	# Full test
	assert result == expected





if __name__ == "__main__":
	pass