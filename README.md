# Supercritical Galton-Watson processes
Let P(z) be a polynomial describing the probability generating function of the offspring distribution of a supercritical Galton-Watson process. The code contained in this repository allows the computation of the density function of the limit distribution W, written as a (truncated) series of generalized Laguerre polynomials multiplied with an exponential function, plus a delta function in 0 corresponding to the probability of extinction.

## Reproducing the experiments from the paper
The script "run_all_tests.m" reproduces all the numerical experiments in the paper [1].

## Computing the limit distribution for a new offspring distribution
The Matlab script "demo.m" shows how to use the function "compute_density_W.m" and how to interpret its outputs.

## References
[1] Alice Cortinovis, Sophie Hautphenne, and Stefano Massei. "Computing the density of the Kesten–Stigum limit in supercritical Galton–Watson processes" (https://arxiv.org/abs/2601.19633).
