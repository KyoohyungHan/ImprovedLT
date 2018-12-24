# Improved HEAAN
The original HEAAN is in https://github.com/snucrypto/HEAAN. This code is for the paper "Faster Homomorphic Discrete Fourier Transforms and Improved FHE Bootstrapping" (https://eprint.iacr.org/2018/1073.pdf).

## Dependency
This library is written by c++ and using NTL library (http://www.shoup.net/ntl/).
We also used Eigen library for sparse matrix operations.

## How to run this code?
First, you need to build HEAAN library in HEAAN/lib folder by the command "make all".
After that, you just need to type "make" and run "test_homomorphic_dft" and "test_improved_bootstrapping" file in each folder in "app".
