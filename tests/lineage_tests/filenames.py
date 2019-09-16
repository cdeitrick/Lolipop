from pathlib import Path
from dataclasses import dataclass
DATA_FOLDER = Path(__file__).parent.parent / "data"
TABLES_FOLDER = DATA_FOLDER / "tables"

generic_tables = [
	#TABLES_FOLDER / "generic.coexistinglineages.xlsx",
	TABLES_FOLDER / "generic.genotypes.3.xlsx",
	TABLES_FOLDER / "generic.genotypes.5.xlsx",
	TABLES_FOLDER / "generic.genotypes.10.xlsx",
	TABLES_FOLDER / "generic.small.xlsx"
]

# Use classes to help with organization.
class ModelTables:
	model_clonal_interferance: Path = TABLES_FOLDER / "model.clonalinterferance.xlsx"
	model_periodic_selection: Path = TABLES_FOLDER / "model.periodicselection.xlsx"
	model_strong_selection: Path = TABLES_FOLDER / "model.strongselection.xlsx"



class RealTables:
	nature_12344 = DATA_FOLDER / "real.nature12344-s2.BYB1-G07.xlsx"


FILENAME_TRUTHSET = DATA_FOLDER /"truthsets" / "truthset.model.area.xlsx"