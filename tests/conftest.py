from pathlib import Path
import pytest
DATA_FOLDER = Path(__file__).parent / "data" / "tables"

generic_tables = [
	DATA_FOLDER / "generic.coexistinglineages.xlsx",
	DATA_FOLDER / "generic.genotypes.3.xlsx",
	DATA_FOLDER / "generic.genotypes.5.xlsx",
	DATA_FOLDER / "generic.genotypes.10.xlsx",
	DATA_FOLDER / "generic.small.xlsx"
]

model_tables = [
	DATA_FOLDER / "model.clonalinterferance.xlsx",
	DATA_FOLDER / "model.periodicselection.xlsx",
	DATA_FOLDER / "model.strongselection.xlsx"
]

real_tables = [
	DATA_FOLDER / "real.nature12344-s2.BYB1-G07.xlsx"
]

@pytest.fixture
def get_generic_tables():
	return generic_tables


def get_model_tables():
	return model_tables


def get_real_tables():
	return real_tables
