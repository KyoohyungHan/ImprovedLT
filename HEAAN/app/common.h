// Standard Headers
#include <algorithm>
#include <math.h>

// Headers for the Eigen lib
#include <Eigen/Sparse>
#include <Eigen/Dense>

// Header for the HEAAN lib
#include "../src/HEAAN.h"

using namespace std;
using namespace Eigen;

typedef SparseMatrix<complex<double>> MatrixXc;

// Conjugation for Sparse Matrix
void conjugate(MatrixXc& A)
{
	for (long k = 0; k < A.outerSize(); ++k) {
		for (MatrixXc::InnerIterator it(A,k); it; ++it) {
			A.coeffRef(it.row(), it.col()) = conj(A.coeffRef(it.row(), it.col()));
		}
	}
	return;
}

// Divide By d for Sparse Matrix
void divideBy(MatrixXc& A, double d)
{
	for (long k = 0; k < A.outerSize(); ++k) {
		for (MatrixXc::InnerIterator it(A,k); it; ++it) {
			A.coeffRef(it.row(), it.col()) = A.coeffRef(it.row(), it.col()) / d;
		}
	}
	return;	
}

// Compare Two array of complex elements with length n
double diff(complex<double>* a1, complex<double>* a2, long n)
{
	double avg = 0.;
	for(long i = 0; i < n; i++) {
		avg += abs(a1[i] - a2[i]);
	}
	avg /= n;
	return log2(avg);
}

uint32_t bitReverse(uint32_t x)
{
	x = (((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1));
	x = (((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2));
	x = (((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4));
	x = (((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8));
	return ((x >> 16) | (x << 16));
}

void bitReverseAndEqual(complex<double>* a, int n)
{
	complex<double>* b = new complex<double>[n];
	for(long i = 0; i < n; i++) {
		long idx = (bitReverse(i) >> (32 - (int)log2(n)));
		b[i] = a[idx];
	}
	for(long i = 0; i < n; i++) {
		a[i] = b[i];
	}
	delete[] b;
}

// this code is copied from "https://en.wikipedia.org/wiki/Cooley–Tukey_FFT_algorithm"
void separate(complex<double>* a, int n)
{
   complex<double>* b = new complex<double>[n / 2];
   for (int i = 0; i<n / 2; i++) { b[i] = a[i * 2 + 1]; }
   for (int i = 0; i<n / 2; i++) { a[i] = a[i * 2]; }
   for (int i = 0; i<n / 2; i++) { a[i + n / 2] = b[i]; }
   delete[] b;
}

// this code is copied from "https://en.wikipedia.org/wiki/Cooley–Tukey_FFT_algorithm"
void FFT(complex<double>* x, int n)
{
   if (n < 2) { }
   else {
      separate(x, n);
      FFT(x, n / 2);
      FFT(x + n / 2, n / 2);
      for (long k = 0; k < n / 2; k++) {
         complex<double> e = x[k];
         complex<double> o = x[k + n / 2];
         complex<double> w = exp(complex<double>(0, -2. * M_PI * k / n));
         x[k] = e + w * o;
         x[k + n / 2] = e - w * o;
      }
   }
}