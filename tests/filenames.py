from pathlib import Path
from loguru import logger
from typing import Dict, Any
DATA_FOLDER = Path(__file__).parent / "data"
TABLES_FOLDER = DATA_FOLDER / "tables"
# These tables have multiple subtables coresponding to different parts of the workflow
# TODO: Add the lang dataset
TRUTHSET_FOLDER = DATA_FOLDER / "truthsets"

generic_tables = {
	#"generic.small": TRUTHSET_FOLDER / "generic.small.xlsx",
	"generic.genotypes.3": TRUTHSET_FOLDER / "generic.genotypes.3.xlsx",
	"generic.genotypes.5": TRUTHSET_FOLDER / "generic.genotypes.5.xlsx",
	"generic.genotypes.10": TRUTHSET_FOLDER / "generic.genotypes.10.xlsx"
}

fake_tables = {
	"generic.genotypes.10": TRUTHSET_FOLDER / "generic.genotypes.10.xlsx",
	# This should only be used to test the area methods.
	"generic.model.area": TRUTHSET_FOLDER / "truthset.model.area.xlsx"
}

# Use classes to help with organization.
model_tables = {
	'model.clonalinterferance': TRUTHSET_FOLDER / "model.clonalinterferance.xlsx",
	'model.periodicselection':  TRUTHSET_FOLDER / "model.periodicselection.xlsx",
	'model.strongselection':   TRUTHSET_FOLDER / "model.strongselection.xlsx"
}

real_tables = {
	'nature12344': TRUTHSET_FOLDER / "real.nature12344-s2.BYB1-G07.xlsx",
	'traverse': TABLES_FOLDER / "real.traverse.tsv",
	'B1': TRUTHSET_FOLDER / "real.B1_muller_try1.xlsx"
}
# These are mainly to test whether Lolipop can load different filetypes.
trajectory_tables = {
	"B1.tsv": DATA_FOLDER / "tables_input_genotypes" / "B1_muller_try1.muller_genotypes.original.tsv",
	"B1.xls": DATA_FOLDER / "tables_input_genotypes" / "B1_muller_try1.muller_genotypes.original.xls",
	"B1.xlsx": DATA_FOLDER / "tables_input_genotypes" / "B1_muller_try1.muller_genotypes.original.xlsx"
}

genotype_tables = {
	"B1.tsv": DATA_FOLDER / "tables_input_genotypes" / "B1_muller_try1.muller_genotypes.original.tsv",
	"B1.xls": DATA_FOLDER / "tables_input_genotypes" / "B1_muller_try1.muller_genotypes.original.xls",
	"B1.xlsx": DATA_FOLDER / "tables_input_genotypes" / "B1_muller_try1.muller_genotypes.original.xlsx"
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

