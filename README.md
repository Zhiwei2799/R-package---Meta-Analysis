# R-package---Meta-Analysis
R package for various meta analysis methods to combine p-values. A meta-analysis is a statistical technique used in research to combine and analyze data from multiple independent studies on a specific topic or research question. Meta-analysis methods typically do not directly combine p-values but instead combine effect sizes or other relevant statistical measures from individual studies. 
## Combine p-values
### Fisher 
The Fisher’s method sums up the logtransformed p-values obtained from individual studies. The combined Fisher’s statistic $X_{fisher}^2 = \sum\limits_{i=1}^{k}\ln(p_i)$ follows a $X^2$ distribution with 2k degrees of freedom under the null hypothesis (assuming null p-values are uniformly distributed).
### Stouffer 
Stouffer's method is used to combine the results from several independence tests bearing upon the same overall hypothesis (H0). The method sums the inverse normal transformed p-values. ${Z_{Stouffer} = \sum\limits_{i=1}^{k} Z_{i}/\sqrt{k}}$
### Maximum p-value (maxP)
The maxP method takes maximum p-value as the test statistic. It follows a beta distribution with degrees of freedom ${\alpha = K}$ and ${\beta = 1}$ under the null hypothesis
### Minimum p-value (minP)
The minP method takes the minimum p-vlaue among the K studies as the test statistic. It follows a beta distribution with degrees of freedom  ${\alpha = 1}$ and ${\beta = k}$ under the null hypothesis.
### Weighted Stouffer
weight is the square root of the sample size, wi
${Z_{weighted~Stouffer}} = \frac{\sum\limits_{i=1}^{k} w_i Z_i}{\sqrt{\sum\limits_{i=1}^{k} w_i^2}}$
### Lancaster
Lancaster's method generalized Fisher's method by assignning different weights using the additivity of chi-squared distribution. ${T = \sum\limits_{i=1}^{k} [X_{n_i}^{2}]^{-1}(1-p_i)}$ where ${[X_{n_i}^{2}]^{-1}}$ is the inverse cumulative chi-square distributon function with ni degrees of freedom

# Reference:  

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-368   

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135688/   
