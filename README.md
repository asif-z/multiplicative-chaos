# multiplicative chaos in number theory

This repository hosts code and data for the paper "A conjectural asymptotic formula for multiplicative chaos in number theory" by D. Aggarwal, U. Subedi, W. Verreault, A. Zaman, and C. Zheng. 

* [generate_A_N_for_standard_complex_Gaussian_samples.cpp](https://github.com/asif-z/multiplicative-chaos/blob/main/generate_A_N_for_standard_complex_Gaussian_samples.cpp) has `C++` code for the case when (X_k) is a sequence of independent standard complex Gaussians.
* [generate_A_N_for_standard_real_Gaussian_samples.cpp](https://github.com/asif-z/multiplicative-chaos/blob/main/generate_A_N_for_standard_real_Gaussian_samples.cpp) has `C++` code for the case when (X_k) is a sequence of independent standard real Gaussians.
* [generate_A_N_for_pm1_samples.cpp](https://github.com/asif-z/multiplicative-chaos/blob/main/generate_A_N_for_pm1_samples.cpp) has `C++` code for the case when (X_k) is a sequence of independent random variables uniform on {\pm 1}.


Our code is written specifically for the Compute Canada clusters, which is where the majority of our computations were performed, and utilizes `OpenMp` and `MPI` for multi-threading.

* [compute_canada.md](https://github.com/asif-z/multiplicative-chaos/blob/main/compute_canada.md) has a brief introduction to using [Compute Canada](https://www.computecanada.ca) clusters.

For example, `generate_A_N_for_standard_complex_Gaussian_samples.cpp` can be compiled by

```
mpiicpc -qopenmp generate_A_N_for_standard_complex_Gaussian_samples.cpp -o complex_Gaussian
```

and can be invoked using 

```
mpirun ./complex_Gaussian
```

This will run a Monte Carlo simulation for (A(n)) in the case when (X_k) is a sequence of independent standard complex Gaussians for N = 20000 and a sample size of 10 million. Intermediate data for A(n) is _not_ saved and only the final sample mean of |A(n)| is written to file.
