#include "new_bootstrapping.h"

struct Parameter
{
	long logN; long logQ;
	long logp; long logc;
	long log2n; long radix;
	long logq;
	long logT;
};

Parameter bootstrapping_test_param1 = {16, 850, 30, 30, 6, 4, 35, 4};
Parameter bootstrapping_test_param2 = {16, 850, 30, 30, 9, 8, 35, 4};
Parameter bootstrapping_test_param3 = {16, 850, 30, 30, 12, 16, 35, 4};

void bootstrapping_test(Parameter parameter)
{
	// HE parameter //
	long logN = parameter.logN;
	long logQ = parameter.logQ;
	long logp = parameter.logp;
	long logc = parameter.logc;

	// Decomposition related parameter //
	long log2n = parameter.log2n;
	long radix = parameter.radix;

	// Bootstrapping parameter //
	long logq = parameter.logq;
	long logT = parameter.logT;

	long n = 1 << log2n;

	cout << "\n***************************" << endl;
	cout << "Test for Improved Bootstrapping" << endl;
	cout << "logN = " << logN << ", logQ = " << logQ << ", logp = " << logp << ", logc = " << logc << endl;
	cout << "slots = " << n << ", radix = " << radix << ", logq = " << logq << ", logT = " << logT << endl;
	cout << "***************************" << endl;
	cout << endl;

	TimeUtils timeutils;
	timeutils.start("KeyGen");
	Ring ring(logN, logQ);
	SecretKey secretKey(ring);
	Scheme scheme(secretKey, ring);
	scheme.addConjKey(secretKey);
	scheme.addLeftRotKeys(secretKey);
	timeutils.stop("KeyGen");

	timeutils.start("Bootstrapping Helper construct");
	BootHelper boothelper(log2n, radix, logc, scheme, ring, secretKey);
	timeutils.stop("Bootstrapping Helper construct");

	complex<double>* mvec = EvaluatorUtils::randomComplexArray(n);
	Ciphertext cipher = scheme.encrypt(mvec, n, logp, logq);
	mvec = scheme.decrypt(secretKey, cipher);

	complex<double>* dvec;
	
	// Bootstrapping //
	timeutils.start("Improved bootstrapping");
	boothelper.bootstrapping(cipher, logq, logQ, logT);
	timeutils.stop("Improved bootstrapping");

	cout << "* Befor logQ = " << parameter.logq << endl;
	cout << "* After logQ = " << cipher.logq << endl;

	// Print Result and Difference //
	dvec = scheme.decrypt(secretKey, cipher);
	cout << "log2(avg of error) = " << diff(mvec, dvec, n) << endl;

	if(mvec != NULL) delete[] mvec;
	if(dvec != NULL) delete[] dvec;
	return;
}

int main(int argc, char* argv[]) {

	bootstrapping_test(bootstrapping_test_param1);
	//bootstrapping_test(bootstrapping_test_param2);
	//bootstrapping_test(bootstrapping_test_param3);

	return 0;
}