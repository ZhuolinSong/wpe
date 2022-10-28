1. The code implements the sequential aggregation method proposed in Zhang and Chen 2020, "Nonparametric covariance estimation for mixed longitudinal studies, with applications in midlife women’s health. " Statistica Sinica.

2. The main functions are in getA.R (detailed descriptions can be found at the beginning of each function). After obtaining the A matrix, the final covariance estimator is $A%*%t(A)$.

3.  The function "getA1_new_eig" uses all, including incomplete, samples to calculate each piece of the A matrix, with or without smoothing.

    The function "getA2_new_eig" uses the possibly incomplete samples that has sufficient coverage to calculate each piece of the A matrix, with/without smoothing.

    The "getA1_new_CV" and "getA2_new_CV" implement the five-fold cross-validation method for the selection of the rank $r$. 