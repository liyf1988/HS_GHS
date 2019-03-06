HS-GHS Version 1.0 03/03/2019

# Description
Draw Monte Carlo samples from the posterior distribution under the graphical horseshoe prior, to estimate the coefficients and the precision matrix in multivariate Gaussian regressions.

# Usage
[beta_save,lambda_sq_save,tau_sq_save,omega_save,lambda_G_sq_save,tau_G_sq_save] = HSGHS(X,Y,1000,5000,eye(q))

beta_est = mean(beta_save,3)

omega_est = mean(omega_save,3)

# Arguments
X: n by p matrix of observed predictors

Y: n by q matrix of observed responses

burning, nmc: number of MCMC burnins and saved samples

Sigma_init: initial value of covariance matrix; need to be symmetric and positive definite q by q matrix

# Value
beta_save: p by q by nmc matrix of saved posterior samples of coefficients

lambda_sq_save: p by q by nmc matrix of saved posterior samples of lambda squared (local tuning parameter), for coefficients

tau_sq_save: 1 by nmc vector of saved posterior samples of tau squared (global tuning parameter), for coefficients

omega_save: q by q by nmc matrix of saved posterior samples of precision matrix

lambda_sq_G_save: q(q-1)/2 by nmc matrix of saved posterior samples of lambda squared (local tuning parameter), for precision matrix

tau_sq_G_save: 1 by nmc matrix of saved posterior samples of tau squared (global tuning parameter), for precision matrix

# Details
Draw posterior samples to estimate the coefficients and the precision matrix in multivariate Gaussian regressions. Posterior mean of the samples is the HS-GHS estimate. The function is a full Gibbs sampler. Within a Gibbs step, the algorithm uses fast sampling method by Bhattacharya et. al. (2016) to sample coefficients, and the graphical horseshoe algorithm by Li et. al. (2017) to sample the precision matrix.

# References
Bhattacharya, Anirban, Antik Chakraborty, and Bani K. Mallick (2016). "Fast sampling with Gaussian scale mixture priors in high-dimensional regression." Biometrika 103, 985-991.

Li, Yunfan, Bruce A. Craig, and Anindya Bhadra (2019). "The Graphical Horseshoe Estimator for Inverse Covariance Matrices." Journal of Computational and Graphial Statistics to appear.

# Examples
See HSGHS_example.m
