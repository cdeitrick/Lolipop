from pathlib import Path
from typing import Dict, Tuple

import pytest

from muller.graphics import Palette, graphicsio


@pytest.fixture
def project_folder(tmp_dir) -> Path:
	return tmp_dir / "project"


@pytest.fixture
def workflow(project_folder) -> graphicsio.OutputGeneratorGraphics:
	palette_1 = None
	palette_2 = None
	palettes = [palette_1, palette_2]

	app = graphicsio.OutputGeneratorGraphics(
		project_folder,
		palettes,
		sample_basename = "simple",
		render = True
	)
	return app


@pytest.fixture
def palettes() -> Tuple[Palette, Palette]:
	genotype_colors_clade: Dict[str,str] = {
		"genotype-1":        "#fca082",
		"genotype-10":       "#ea0000",
		"genotype-11":       "#abd0e6",
		"genotype-2":        "#c6c6c6",
		"genotype-3":        "#e32f27",
		"genotype-4":        "#3787c0",
		"genotype-5":        "#686868",
		"genotype-6":        "#796eb2",
		"genotype-7":        "#ffca00",
		"genotype-8":        "#73c476",
		"genotype-9":        "#c6c7e1",
		"genotype-filtered": "#333333",
		"genotype-0":        "#FFFFFF"
	}
	genotype_colors_unique: Dict[str,str] = {
		"genotype-1":        "#e6194b",
		"genotype-10":       "#fabebe",
		"genotype-11":       "#008080",
		"genotype-2":        "#3cb44b",
		"genotype-3":        "#ffe119",
		"genotype-4":        "#4363d8",
		"genotype-5":        "#f58231",
		"genotype-6":        "#911eb4",
		"genotype-7":        "#46f0f0",
		"genotype-8":        "#f032e6",
		"genotype-9":        "#bcf60c",
		"genotype-filtered": "#333333",
		"genotype-0":        "#FFFFFF"
	}

	p1 = Palette("clade", genotype_colors_clade)
	p2 = Palette("unique", genotype_colors_unique)
	return p1, p2

def test_outputgraphicsworkflow_with_various_options(tmp_path, palettes):
	""" basically tests if__init__() craches."""

	folder = tmp_path / "project"


if __name__ == "__main__":
	pass
