from pathlib import Path

import matplotlib.pyplot as plt
import pandas


# plt.switch_backend('QT4Agg')


class Plot:
	def __init__(self):
		self.xlabel = "Time"
		self.ylabel = "Frequency"
		self.axlabelsize = 24

		self.linewidth = 10
		self.axticklabelsize = '20'

	def _formatplot(self, ax: plt.Axes) -> plt.Axes:
		ax.set_xlabel(self.xlabel, fontsize = self.axlabelsize)
		ax.set_ylabel(self.ylabel, fontsize = self.axlabelsize)

		# Increase the font size of the frequency and timepoint labels
		for label in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
			label.set_fontsize(self.axticklabelsize)

		ax.legend()
		return ax

	@staticmethod
	def _extract_color(label: str) -> str:
		elements = label.split('-')
		color = elements[1]  # Should always be the second element.
		return color

	def run(self, table: pandas.DataFrame, filename: Path):
		if 'Genotype' in table.columns:  # Should be the index for genotype tables
			del table['Genotype']

		fig, ax = plt.subplots(figsize = (12, 10))

		seen_labels = set()
		for index, row in table.iterrows():
			color = self._extract_color(index)
			genotype_label = '-'.join(index.split('-')[:2])
			label = genotype_label if genotype_label not in seen_labels else None
			seen_labels.add(genotype_label)
			ax.plot(row.index, row.values, label = label, color = color, linewidth = self.linewidth)

		ax = self._formatplot(ax)
		ax.set_xlim(0, max(row.index))
		ax.set_ylim(0, 1)
		plt.tight_layout()
		plt.savefig(str(filename), dpi = 250)


if __name__ == "__main__":
	# Should have both genotypes and trajectories defined.

	folder = Path(__file__).parent / "tests" / "data"
	tables_folder = folder / "tables"
	table_filename = tables_folder / "model.clonalinterferance.xlsx"

	trajectories = pandas.read_excel(table_filename, sheet_name = "trajectory").set_index('Trajectory')
	genotypes = pandas.read_excel(table_filename, sheet_name = "genotype").set_index('Genotype')

	plotter = Plot()
	plotter.run(trajectories, table_filename.with_suffix(f'.trajectory.png'))
	plotter.run(genotypes, table_filename.with_suffix(".genotype.png"))
