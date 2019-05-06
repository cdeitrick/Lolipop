 # Nesting successive genotypes

## Checks
After grouping mutational trajectories into genotypes, the scripts then attempt which genotypes arise in the background of which other genotypes.
Genotypes are nested according to a few basic rules:
### Summation
Check whether the unnested and nested genotype consistently sum to greater than 1. This is weak evidence that one of the genotypes is a background and the other arises in the background of the first. The implementation currently checks whether the sum of both trajectories is greater than 1.15 (chosen based on the original implementation) at least once, or greater than 1.03 (based on the user-specified uncertainty option) at least twice. Future versions of the script will implement this consistency check in a more statistically-relevant way.

### Size
This checks Whether one of the genotypes is consistently larger than the other. A genotype is considered consistently greater than another if the difference between the two frequencies exceeds a value of 0.15 (based on original scripts) at least once or exceeds 0.03 (based on the uncertainty option) at least twice. Future versions of the script will define this in a more statistically-relevant manner.

### Correlation
The covariance between two series of random variables provides a measure of the correlation between both series. This is also used to calculate the pearson distance metric. A genotype which arises in the background of another genotype must be correlated with the parent genotype.
The covariance of two series $X$ and $Y$ is defined as
$$
Cov(X,Y) = \frac{1}{n}\sum_{i=0}^n (X_i-\bar{X})(Y_i-\bar{Y})
$$

### The Jaccard Distance

The jaccard distance is a measurement of the dissimilarity of two sample sets. Since each frequency measurement describes the abundance of a mutation at a sampled timepoint, we can use the jaccard distance to compare the abundance of two genotypes over the course of an experiment.
Since each genotype represents a set of abundance measurements, the jaccard distance between two genotypes $X$ and $Y$ can be calculated as follows:
$$
J(X,Y)=1 - \frac{|X\cup Y|}{|X|+|Y|+|X \cap Y|}\equiv \frac{|X\cup Y|-|X \cap Y|}{|X \cup Y|}
$$
Since each genotype is a set of frequency measurements, the cardinality of each set can be computed as
$$
|X| = \sum_{i=0}^n X_i
$$


Based on this, the union and intersection can be defined as
$$
|X \cap Y| = \sum_{i=0}^n min(X_i,Y_i)
$$
$$
|X \cup Y| = |X| + |Y| - |X \cap Y|
$$


When $Y$ arises in $X$, all elements of $Y$ are also contained in $X$, $|X \cup Y| = |X|$ and $|X \cap Y| = |Y|$, reducing the above equation to
$$
J(X,Y) = \frac{|X|-|Y|}{|X|}
$$


If the computed jaccard distance is equal to the above equation, genotype $Y$ arises in genotype $X$.

#### Example

| Genotype   | 0 | 17 | 25 | 44    | 66    | 75    | 90    |
|------------|---|----|----|-------|-------|-------|-------|
| genotype-6 | 0 | 0  | 0  | 0.273 | 0.781 | 1.000 | 1.000 |
| genotype-7 | 0 | 0  | 0  | 0.403 | 0.489 | 0.057 | 0.080 |
| genotype-3 | 0 | 0  | 0  | 0     | 0.211 | 0.811 | 0.813 |

The jaccard distance between genotype-6 and genotype-7, assuming genotype-6 is the background:
$$|X| = 0.273+0.781+1+1 = 3.054$$
$$|Y| = 0.403+0.489+0.057+0.080=1.029$$
$$|X \cap Y|=.273+.489+.057+0.080=0.899$$
$$|X \cup Y| \equiv |X|+|Y|-|X \cap Y|= 3.054+1.029-0.899=3.184$$
Now let's calculate the jaccard distance the traditional way:
$$
J(X,Y)=\frac{|X\cup Y|-|X \cap Y|}{|X \cup Y|}=\frac{3.184-0.899}{3.184}=0.718
$$
Test against the ideal distance:
$$
J(X,Y) = \frac{|X|-|Y|}{|X|}=\frac{3.054-1.029}{3.054}=0.663
$$

Since $0.718\ne 0.663$, We cannot say that $Y$ arises in the background of $X$.

Let's test if genotype-3 arises in the background of genotype-6:

$$|X| = 0.273+0.781+1.000+1.000 = 3.054$$
$$|Y| = 0.211+0.811+0.813=1.835$$
$$|X \cap Y|=0+.211+.811+.813=1.835$$
$$|X \cup Y| \equiv |X|+|Y|-|X \cap Y|= 3.054+1.835-1.835=3.184$$

The jaccard distance is then
$$
J(X,Y)=\frac{|X\cup Y|-|X \cap Y|}{|X \cup Y|}=\frac{3.184-1.835}{3.184}=0.425
$$
Test against the ideal distance:
$$
J(X,Y) = \frac{|X|-|Y|}{|X|}=\frac{3.054-1.835}{3.054}=0.425
$$

So, there is evidence that genotype-3 arises in the background of genotype-6.



## Ruleset

Given the above checks, this is the ruleset for determining if a nested genotype is a potential background for the unnested genotype being tested.

1. If the genotypes are negatively correlated, the nested genotype is not a background candidate for the current genotype.
2. If the unnested genotype contains more elements than the nested genotype, it cannot be a background candidate.
3. If the nested genotype is consistently larger than the unnested genotype, is positively correlated with the unnested genotype, and the jaccard distance indicates the unnested genotype is a subset of the nested genotype, then the nested genotype is assigned as the background for the unnested genotype. The newly nested genotype is then available as a potential background candidate for subsequent unnested genotypes.
4. If the unnested genotype is consistently smaller than the nested genotype or is determined to be a subset of the nested genotype, the nested genotype is considered a potential background for the unnested genotype. The unnested genotype will continue to be tested against all other nested genotypes. If no other genotype has greater evidence for being the background of the unnested genotype, the script prioritises the nested genotype with the greatest maximum frequency which was detected closest to when the unnested genotype was first detected.
