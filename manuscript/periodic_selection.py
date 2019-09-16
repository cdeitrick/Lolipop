from pathlib import Path 
import pandas
import math
def cdf(x, mean: float = 0.5, variance: float = 0.1, offset = 0):
	x += offset
	A = (x - mean) / (variance * math.sqrt(2))
	result = 0.5 * (1 + math.erf(A))

	return result

if __name__ == "__main__":
	import numpy 
	x = numpy.linspace(-4,4,100)
	y = [0.5*(1+math.erf(i) for i in x)]

	import matplotlib.pyplot as plt 
	plt.plot(x, y)
	plt.show()