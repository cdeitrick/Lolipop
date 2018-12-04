

NOTE: This is incorrect. Kept for reference.

The script calculates the relative similarity between all pairs of mutations with frequencies $0<=f_i<=1$. These are consistent with $n$ independent draws from a normal distribution, assuming an average variance of  

1. $$\bar{f_i} = (f_{ai}+f_{bi})/2$$

2. $$\sigma_i^2 = n_{binom}\bar{f_i}(1-\bar{f_i})$$

where $\bar{f_i}$ is the average of the two frequencies at each time point, and $n_{binom}$ is picked arbitrarily as a value that gives reasonable uncertainties given our data (currently $1/5$). For a time series with $n_p$ random draws (observations) will have

3. $\sigma_{tot}^2 = \sum_i^{n_p} \sigma_{i}^2 = n_{binom} \sum_i^{n_p} \bar{f_i}(1 - \bar{f_i}) $

from a property of normal distributions.

For the mean of all frequencies in the averaged time series, 

4. $\bar{X} = \sum_{i}^{n_X} X_i/n_X = 1/n_X \sum_{i}^{n_X} X_i$
5. $\bar{f} = \sum_{i}^{n_p} f_i/n_p$

where $X$ is the mean frequency of both time series at each time point and $n_X$ is the number of time points,

$$\sigma_{\bar{X}} = \sigma_{tot}/n_X = \sqrt( \sum_i f_i(1-f_i))/\sqrt{n_X}$$

Finally, given $\sigma$ and $\bar{X}$ we can construct a p-value for
the measurement by numerically computing the following integral:

$$
1 - \int_{- \bar{X}/ \sigma_{\bar{X}}}^{\bar{X}/\sigma_{\bar{X}}} \frac{1}{\sqrt{2 \pi}} e^{-x^2/2} dx
$$
Which is based on the [error function](https://en.wikipedia.org/wiki/Error_function):
$$
erf(x) = \frac{1}{\sqrt{\pi}} \int_{-x}^{x}e^{-t^2}dt
$$

these are consistent with n independent draws from a normal distribution,
assuming an average variance of  $n_{binom} p (1-p)$, where p is the
average of the two frequencies, and $n_{binom}$ is picked arbitrarily as a
value that gives reasonable uncertainties given our data

n random draws will have
$\sigma_{tot}^2 = \sum_i \sigma_i^2 =  n_{binom} \sum_i p_i(1 - p_i)$, from
a property of normal distributions

for $\bar{X} = \sum_{i}^{n_X} X/n_X$

$\sigma_{\bar{X}} = \sigma_{tot}/n_X = \sqrt( \sum_i p_i(1-p_i))/\sqrt{n_X}$

Finally, given \sigma and \bar{X} we can construct a p-value for
the measurement by numerically computing the following integral:

$$ 1 - \int_{- \bar{X}/ \sigma_{\bar{X}}}^{\bar{X}/ \sigma_{\bar{X}}} \frac{1}{\sqrt{2 \pi}} e^{-x^2/2} dx$$

# Code 

$p_s = (p_{ai} + p_{bi})/2$

$\sigma_i^2 = p(1-p)/5$

$d = (p_{ai}-p_{bi})$

$\sigma_{pair} = \sqrt{\sum \sigma_{i}^2} / n_d$

$\bar{d} = \sum d / n_d$

$p = 1-erf(\bar{d}/\sqrt{2}\sigma_{pair})$

$\bar{d}/\sqrt{2}\sigma_{pair} = \sum{d/n_d}/\sqrt{2}\sqrt{\sum \sigma_{i}^i} / n_d$

