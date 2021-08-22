# random-power-series-coefficients
An algorithm for computing random exponential power series coefficients

This repository hosts code and data for the paper "A conjectural asymptotic formula for multiplicative chaos in number theory" by D. Aggarwal, U. Subedi, W. Verreault, A. Zaman, and C. Zheng. 

* [generate_A_N_for_standard_complex_Gaussian_samples.cpp](https://github.com/asif-z/multiplicative-chaos/blob/main/generate_A_N_for_standard_complex_Gaussian_samples.cpp) has C++ code for our algorithm in the case when (X_k) is a sequence of independent standard complex Gaussians.
* [generate_A_N_for_standard_real_Gaussian_samples.cpp](https://github.com/asif-z/multiplicative-chaos/blob/main/generate_A_N_for_standard_real_Gaussian_samples.cpp) has C++ code for our algorithm in the case when (X_k) is a sequence of independent standard real Gaussians.
* [generate_A_N_for_pm1_samples.cpp](https://github.com/asif-z/multiplicative-chaos/blob/main/generate_A_N_for_pm1_samples.cpp) has C++ code for our algorithm in the case when (X_k) is a sequence of independent random variables uniform on {\pm 1}.
* [compute_canada.md](https://github.com/asif-z/multiplicative-chaos/blob/main/compute_canada.md) has a short introduction to using Compute Canada clusters, where the majority of our computations were performed.
