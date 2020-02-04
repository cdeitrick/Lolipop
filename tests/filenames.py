from pathlib import Path
from loguru import logger
from typing import Dict, Any
DATA_FOLDER = Path(__file__).parent / "data"
TABLES_FOLDER = DATA_FOLDER / "tables"

generic_tables_with_trajectories = {
	"generic.genotypes.3": TABLES_FOLDER / "generic.genotypes.3.xlsx",
	"generic.genotypes.5": TABLES_FOLDER / "generic.genotypes.5.xlsx",
	"generic.genotypes.10": TABLES_FOLDER / "generic.genotypes.10.xlsx"
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

# Add a method to compare dictionaries, since pytest's default view makes finding the mismatched key-value a hassle.
def _compare_dictionaries_exclusive_keys(left:Dict[Any, Any], right:Dict[Any,Any]):
	# Find the keys that are exclusive to each dictionary
	exclusive_left = set(left.keys()) - set(right.keys())
	exclusive_right = set(right.keys()) - set(left.keys())

	if len(exclusive_left) > 0:
		logger.debug(f"The folloowing keys are only present in the left dictionary:")
		for el in exclusive_left:
			logger.debug(f"\t{el} -> {left[el]}")
	if len(exclusive_right) > 0:
		logger.debug(f"The following keys are only present in the right dictionary:")
		for er in exclusive_right:
			logger.debug(f"\t{er} -> {left[er]}")
def _compare_dictionaries_conflicting_values(left:Dict[Any, Any], right:Dict[Any,Any]):

	allkeys = set(left.keys()) | set(right.keys())
	for key in sorted(allkeys):
		vl = left.get(key)
		vr = right.get(key)

		if vl != vr:
			logger.debug(f"Found missmatched values for key '{key}': recieved '{vl}' but expected '{vr}'")

def compare_dictionaries(left:Dict[Any, Any], right:Dict[Any,Any]):

	# Find the keys that are exclusive to each dictionary
	_compare_dictionaries_exclusive_keys(left, right)

	# Fins keys with unequal values
	_compare_dictionaries_conflicting_values(left, right)

