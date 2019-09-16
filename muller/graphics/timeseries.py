import matplotlib.pyplot as plt
import pandas
from muller import widgets

class DerivativePlot:
	def __init__(self):
		pass
	def run(self, left:pandas.Series, right: pandas.Series):
		detected_left, detected_right = widgets.get_valid_points(left, right, 0.03, flimit = 0.97, inner = True)

		diff_left = detected_left.diff()
		diff_right = detected_right.diff()

		fig, ax = plt.subplots(figsize = (15,15))
		ax.scatter(diff_left.values, diff_right.values)

		ax.set_xlabel("left derivative")
		ax.set_ylabel("right derivative")

		ax.set_xlims(-1, 1)
		ax.set_ylims(-1, 1)

		plt.show()

if __name__ == "__main__":
	pass


