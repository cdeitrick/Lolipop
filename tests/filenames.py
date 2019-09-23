from pathlib import Path

DATA_FOLDER = Path(__file__).parent / "data"
TABLES_FOLDER = DATA_FOLDER / "tables"

generic_tables_with_trajectories = {
	"generic.genotypes.3": TABLES_FOLDER / "generic.genotypes.3.xlsx",
	"generic.genotypes.5": TABLES_FOLDER / "generic.genotypes.5.xlsx",
}

fake_tables = {
	"generic.genotypes.10": TABLES_FOLDER / "generic.genotypes.10.xlsx",
}

# Use classes to help with organization.
model_tables = {
	'model.clonalinterferance': TABLES_FOLDER / "model.clonalinterferance.xlsx",
	'model.periodicselection':  TABLES_FOLDER / "model.periodicselection.xlsx",
	'model.strongselection':   TABLES_FOLDER / "model.strongselection.xlsx"
}

real_tables = {
	'nature12344': TABLES_FOLDER / "real.nature12344-s2.BYB1-G07.xlsx"
}

FILENAME_TRUTHSET = DATA_FOLDER / "truthsets" / "truthset.model.area.xlsx"
