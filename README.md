# Estimate joint density of &mu; and &sigma^2;
---
Store the "density_estimation" folder in your home folder and run it. In each iteration, we ran 40 simulations in parallel.

**DPMMfunction.R**: R function for estimating mixture of normal-inverse gamma and all other variantions.
(Our Algorithm)

**SUREmethods.R**: R function for estimating means for several other SURE algoritms.
(SURE Estimates for a Heteroscedastic Hierarchical
Model-Xianchao Xie , S. C. Kou & Lawrence D. Brown, On SURE-Type Double Shrinkage Estimation-
Bing-Yi Jing, Zhouping Li, Guangming Pan & Wang Zhou)

**grouplinearfunction_all.R**: R function for estimating means for group linear algorithms.
(Group-Linear Empirical Bayes Estimates for a Heteroscedastic
Normal Mean-Asaf Weinstein, Zhuang Ma, Lawrence D. Brown, Cun-Hui Zhang)

**bash_script.sh**: Bash script to execute all codes below.

**diffq_exam1.R - diffq_exam6.R**: 6 example
q=100,200,...,1000, B<sub>q</sub>=100, X<sub>ij</sub> = &mu;<sub>i</sub> + &sigma;<sub>ij</sub> &epsilon;<sub>i</sub>, &sigma;<sub>i</sub> unknown. We are estimating the joint density of &mu; and &sigma^2;. 