#include "new_homDFT.h"

struct Parameter
{
	long logN; long logQ;
	long logp; long logc;
	long log2n; long radix;
};

Parameter homomorphic_dft_test_param1 = {14, 125, 30, 25, 12, 16};
Parameter homomorphic_dft_test_param2 = {14, 150, 30, 25, 12, 8};
Parameter homomorphic_dft_test_param3 = {14, 210, 30, 25, 12, 4};
Parameter homomorphic_dft_test_param4 = {15, 350, 30, 25, 12, 2};

void homomorphic_dft_test(Parameter parameter)
{

	cout << "\n********************************" << endl;
	cout << "********************************" << endl;
	cout << "Test for Homomorphic DFT" << endl;
	cout << "logN = " << parameter.logN << ", logQ = " << parameter.logQ << ", logp = " << parameter.logp << ", logc = " << parameter.logc << endl;
	cout << "log2(size) = " << parameter.log2n << ", radix = " << parameter.radix << endl;
	cout << "********************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(parameter.logN, parameter.logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("HomDFT construct");
	DFTHelper dfthelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
	timeutils.stop("HomDFT construct");

	long n = 1 << parameter.log2n;

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);

	complex<double>* dvec;
	complex<double>* dftvec;
	
	// Homomorphic DFT //
	timeutils.start("homomorphic DFT_NR");
	dfthelper.DFT_NR(cipher);
	timeutils.stop("homomorphic DFT_NR");

	// Plaintext FFT //
	dftvec = new complex<double>[n]();
	for(long i = 0; i < n; i++) dftvec[i] = mvec[i];
	FFT(dftvec, n);	
	bitReverseAndEqual(dftvec, n);
	
	// Print Result and Difference //
	dvec = scheme.decrypt(secretKey, cipher);
	cout << "log2(avg of error) = " << diff(dftvec, dvec, n) << endl;

	if(mvec != NULL) delete[] mvec;
	if(dftvec != NULL) delete[] dftvec;
	if(dvec != NULL) delete[] dvec;
	return;
}

void homomorphic_inv_dft_test(Parameter parameter)
{

	cout << "\n********************************" << endl;
	cout << "********************************" << endl;
	cout << "Test for Homomorphic inv DFT" << endl;
	cout << "logN = " << parameter.logN << ", logQ = " << parameter.logQ << ", logp = " << parameter.logp << ", logc = " << parameter.logc << endl;
	cout << "log2(size) = " << parameter.log2n << ", radix = " << parameter.radix << endl;
	cout << "********************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(parameter.logN, parameter.logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("HomDFT construct");
	DFTHelper dfthelper(parameter.log2n, parameter.radix, parameter.logc, scheme, ring, secretKey);
	timeutils.stop("HomDFT construct");

	long n = 1 << parameter.log2n;

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher = scheme.encrypt(mvec, n, parameter.logp, parameter.logQ);

	complex<double>* dvec;
	complex<double>* dftvec;
	
	// Homomorphic inv DFT //
	timeutils.start("homomorphic inv DFT_NR");
	dfthelper.IDFT_RN(cipher);
	timeutils.stop("homomorphic inv DFT_NR");

	// Plaintext inv FFT //
	dftvec = new complex<double>[n]();
	for(long i = 0; i < n; i++) dftvec[i] = mvec[i];
	bitReverseAndEqual(dftvec, n);
	for(long i = 0; i < n; i++) dftvec[i] = conj(dftvec[i]);
	FFT(dftvec, n);	
	for(long i = 0; i < n; i++) dftvec[i] = conj(dftvec[i]);
	for(long i = 0; i < n; i++) {
		dftvec[i].real(dftvec[i].real() / (double(n)));
		dftvec[i].imag(dftvec[i].imag() / (double(n)));
	}
	
	// Print Result and Difference //
	dvec = scheme.decrypt(secretKey, cipher);
	cout << "log2(avg of error) = " << diff(dftvec, dvec, n) << endl;

	if(mvec != NULL) delete[] mvec;
	if(dftvec != NULL) delete[] dftvec;
	if(dvec != NULL) delete[] dvec;
	return;
}

int main(int argc, char* argv[]) {

	homomorphic_dft_test(homomorphic_dft_test_param1);
	homomorphic_dft_test(homomorphic_dft_test_param2);
	homomorphic_dft_test(homomorphic_dft_test_param3);
	homomorphic_dft_test(homomorphic_dft_test_param4);

	homomorphic_inv_dft_test(homomorphic_dft_test_param1);
	homomorphic_inv_dft_test(homomorphic_dft_test_param2);
	homomorphic_inv_dft_test(homomorphic_dft_test_param3);
	homomorphic_inv_dft_test(homomorphic_dft_test_param4);

	return 0;
}