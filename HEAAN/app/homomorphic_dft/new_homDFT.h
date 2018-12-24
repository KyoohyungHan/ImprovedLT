#include "../common.h"

class DFTHelper
{

public:

	// ring class for encoding
	Ring& ring;

	// scheme class for homomorphic computation
	Scheme& scheme;

	// parameters
	long n;
	long radix;
	long log2n;
	long logrn;
	long log2r;
	long logp;

	// decomposed FFT matrix
	MatrixXc* DFTDecomp;
	MatrixXc* invDFTDecomp;

	// pre-encoded polys
	ZZ*** encodeDFTpoly;
	ZZ*** encodeinvDFTpoly;
	
	DFTHelper(long log2n, long radix, long logp, Scheme& scheme, Ring& ring, SecretKey& secretKey);

	~DFTHelper();

	void buildPartialMatrix(MatrixXc& A, long start_row, long start_col, long size);

	// Build DFT decomposion matrix using Eigen library
	void buildDFTDecompMatrix(MatrixXc* DFTDecomp, long n);

	// Encode diagonal vector of the given matrix A
	void encodeDiagonal(ZZ* poly, MatrixXc A, long idx);

	// Homomorphic DFT (normal order -> reversed order)
	void DFT_NR(Ciphertext& cipher);

	// Homomorphic inverse DFT (reversed order -> normal order)
	void IDFT_RN(Ciphertext& cipher);

};

DFTHelper::DFTHelper(long log2n, long radix, long logp, Scheme& scheme, Ring& ring, SecretKey& secretKey)
	 : scheme(scheme), ring(ring), log2n(log2n), radix(radix), logp(logp)
{
	// parameters		
	n = 1 << log2n;
	log2r = log2(radix);
	logrn = log2n / log2r;
	
	// decompose DFT matrix and invDFT
	DFTDecomp = new MatrixXc[log2n]();
	invDFTDecomp = new MatrixXc[log2n]();
	buildDFTDecompMatrix(DFTDecomp, n);
	for(long i = 0; i < log2n; i++) {
		invDFTDecomp[i] = DFTDecomp[i].transpose();
		conjugate(invDFTDecomp[i]);
		divideBy(invDFTDecomp[i], 2.);
	}
	
	// make public keys for left rotations
	for(long i = 0; i < logrn; i++) {
		long gap = n / (radix * pow(radix, i));
		if(i == 0) {
			if(radix > 4) {
				long bs = 1 << (long)floor((log2(radix)) / 2.);
				long gs = (radix) / bs;
				for(long j = 1; j < bs; j++) {
					long idx = j * gap;
					idx %= n;
					if(idx != 0) scheme.addLeftRotKey(secretKey, idx);
				}
				for(long j = 1; j < gs; j++) {
					long idx = j * gap * bs;
					idx %= n;
					if(idx != 0) scheme.addLeftRotKey(secretKey, idx);
				}
			}
			else {
				for(long j = 1; j < radix; j++) {
					long idx = j * gap;
					idx %= n;
					if(idx != 0) scheme.addLeftRotKey(secretKey, idx);
				}
			}
		}
		else {
			scheme.addRightRotKey(secretKey, (radix - 1) * gap);
			if(radix > 2) {
				long bs = 1 << (long)floor((log2(radix) + 1) / 2.);
				long gs = (2 * radix) / bs;
				for(long j = 1; j < bs; j++) {
					long idx = j * gap;
					idx %= n;
					if(idx != 0) scheme.addLeftRotKey(secretKey, idx);
				}
				for(long j = 1; j < gs; j++) {
					long idx = j * gap * bs;
					idx %= n;
					if(idx != 0) scheme.addLeftRotKey(secretKey, idx);
				}
			}
			else {
				for(long j = 1; j < 2 * radix - 1; j++) {
					long idx = j * gap;
					idx %= n;
					if(idx != 0) scheme.addLeftRotKey(secretKey, idx);
				}
			}
		}
	}
	
	// encode diagonals for DFT //
	encodeDFTpoly = new ZZ**[logrn];
	for(long i = 0; i < logrn; i++) {
		for(long j = 1; j < log2r; j++)
			DFTDecomp[i * log2r] = DFTDecomp[i * log2r + j] * DFTDecomp[i * log2r];
		if(i == 0) {
			encodeDFTpoly[0] = new ZZ*[radix];
			for(long j = 0; j < radix; j++) {
				encodeDFTpoly[0][j] = new ZZ[ring.N]();
				long gap = n / radix;
				encodeDiagonal(encodeDFTpoly[0][j], DFTDecomp[0], j * gap);
			}
		} else {
			encodeDFTpoly[i] = new ZZ*[2 * radix - 1];
			long gap = n / (radix * pow(radix, i));
			for(long j = 0; j < 2 * radix - 1; j++) {
				encodeDFTpoly[i][j] = new ZZ[ring.N]();
				encodeDiagonal(encodeDFTpoly[i][j], DFTDecomp[i * log2r], (-radix + 1 + j) * gap);
			}
		}
	}

	// encode diagonals for invDFT //
	encodeinvDFTpoly = new ZZ**[logrn];
	for(long i = 0; i < logrn; i++) {
		MatrixXc tmpMat = DFTDecomp[i * log2r].transpose();
		conjugate(tmpMat);
		divideBy(tmpMat, (double)radix);
		if(i == 0) {
			encodeinvDFTpoly[0] = new ZZ*[radix];
			for(long j = 0; j < radix; j++) {
				encodeinvDFTpoly[0][j] = new ZZ[ring.N]();
				long gap = n / radix;
				encodeDiagonal(encodeinvDFTpoly[0][j], tmpMat, j * gap);
			}
		} else {
			encodeinvDFTpoly[i] = new ZZ*[2 * radix - 1];
			long gap = n / (radix * pow(radix, i));
			for(long j = 0; j < 2 * radix - 1; j++) {
				encodeinvDFTpoly[i][j] = new ZZ[ring.N]();
				encodeDiagonal(encodeinvDFTpoly[i][j], tmpMat, (-radix + 1 + j) * gap);
			}
		}
	}

	// pre-compute
	ZZ* tmpPoly = new ZZ[ring.N](); 
	for(long i = 0; i < logrn; i++) {
		long bs, gs;
		long gap = n / (radix * pow(radix, i)); 
		if(i == 0) {
			if(radix > 4) {
				long bs = 1 << (long)floor((log2(radix)) / 2.);
				long gs = (radix) / bs;
				for(long j = 1; j < gs; ++j) {
					for(long k = 0; k < bs; ++k) {
						ring.leftRotate(tmpPoly, encodeDFTpoly[0][j*bs+k],n-j*bs*gap);
						for(long l = 0; l < ring.N; ++l)
							encodeDFTpoly[0][j*bs+k][l] = tmpPoly[l];
						ring.leftRotate(tmpPoly, encodeinvDFTpoly[0][j*bs+k],n-j*bs*gap);
						for(long l = 0; l < ring.N; ++l)
							encodeinvDFTpoly[0][j*bs+k][l] = tmpPoly[l];
					}
				}
			}
		} 
		else {
			if(radix > 2) {
				long bs = 1 << (long)floor((log2(radix) + 1) / 2.);
				long gs = (2 * radix) / bs;
				for(long j = 1; j < gs; ++j) {
					for(long k = 0; k < bs; ++k) {
						if((j * bs + k) < (2 * radix - 1)) {
							ring.leftRotate(tmpPoly, encodeDFTpoly[i][j*bs+k],n-j*bs*gap);
							for(long l = 0; l < ring.N; ++l)
								encodeDFTpoly[i][j*bs+k][l] = tmpPoly[l];
							ring.leftRotate(tmpPoly, encodeinvDFTpoly[i][j*bs+k],n-j*bs*gap);
							for(long l = 0; l < ring.N; ++l)
								encodeinvDFTpoly[i][j*bs+k][l] = tmpPoly[l];
						}
					}
				}
			}
		}
	}
	delete[] tmpPoly;
}

void DFTHelper::buildPartialMatrix(MatrixXc& A, long start_row, long start_col, long size)
{
	long hsize = size / 2;
	complex<double> omega = exp(complex<double>(0, -2 * M_PI / size));
	for(long i = 0; i < hsize; i++) {
		A.coeffRef(start_row + i, start_col + i) = 1.;
		A.coeffRef(start_row + hsize + i, start_col + i) = pow(omega, i);;
		A.coeffRef(start_row + i, start_col + hsize + i) = 1.;
		A.coeffRef(start_row + hsize + i, start_col + hsize + i) = -pow(omega, i);
	}
	return;
}

void DFTHelper::buildDFTDecompMatrix(MatrixXc* DFTDecomp, long n)
{
	long log2n = (long)log2(n);
	for(long i = 0; i < log2n; i++) {
		DFTDecomp[i].resize(n, n);
		for(long j = 0; j < (1 << i); j++) {
			long size = n / (1 << i);
			buildPartialMatrix(DFTDecomp[i], size * j, size * j, size);
		}
	}
	return;
}

void DFTHelper::encodeDiagonal(ZZ* poly, MatrixXc A, long idx)
{
	long col = A.cols();
	long row = A.rows();
	complex<double>* vec = new complex<double>[row]();
	for(long i = 0; i < row; i++) {
		long tmp = i + idx;
		while(1) {
			if(tmp < 0) tmp += col;
			else break;
		}
		tmp %= col;
		vec[i] = A.coeff(i, tmp);
	}
	ring.encode(poly, vec, row, logp);
	delete[] vec;
	return;
}

void DFTHelper::DFT_NR(Ciphertext& cipher)
{
	// temperary things //
	ZZ* tmpPoly = new ZZ[ring.N]();
	Ciphertext tmpCipher;

	// first step of DFT //
	// baby-step & giant-step method //
	long gap = n / radix;
	long bs = 1 << (long)floor((log2(radix)) / 2.);
	long gs = (radix) / bs;
	if(radix > 4) {
		// In this case, we use basy-step giant-step method
		Ciphertext* bsrotFirst = new Ciphertext[bs];
		Ciphertext* gsrotFirst = new Ciphertext[gs];		
		bsrotFirst[0] = cipher;
		for(long i = 0; i < bs-1; ++i) { // baby-step
			bsrotFirst[i+1] = scheme.leftRotateFast(cipher, (i+1)*gap);
		}
		for(long i = 0; i < gs; i++) {
			gsrotFirst[i] = scheme.multByPoly(bsrotFirst[0], encodeDFTpoly[0][i*bs], logp);
			for(long j = 1; j < bs; j++) {
				Ciphertext tmp = scheme.multByPoly(bsrotFirst[j], encodeDFTpoly[0][i*bs+j], logp);
				scheme.addAndEqual(gsrotFirst[i], tmp);
			}
		}
		cipher = gsrotFirst[0];
		for(long i = 0; i < gs-1; ++i) { // giant-step
			scheme.leftRotateFastAndEqual(gsrotFirst[i+1], (i+1)*bs*gap);
		}
		for(long i = 1; i < gs; i++) {
			scheme.addAndEqual(cipher, gsrotFirst[i]);
		}
		delete[] bsrotFirst;
		delete[] gsrotFirst;
	} else {
		tmpCipher = cipher;
		scheme.multByPolyAndEqual(cipher, encodeDFTpoly[0][0], logp);
		for(long i = 1; i < radix; i++) {
			long gap = n / radix;
			Ciphertext rot = scheme.leftRotateFast(tmpCipher, i * gap);
			scheme.multByPolyAndEqual(rot, encodeDFTpoly[0][i], logp);
			scheme.addAndEqual(cipher, rot);
		}
	}
	scheme.reScaleByAndEqual(cipher, logp);

	// other step for FFT //
	if(radix > 2) {
		// In this case, we use baby-step giant-step method //
		bs = 1 << (long)floor((log2(radix) + 1) / 2.);
		gs = (2 * radix) / bs;
		Ciphertext* bsrot = new Ciphertext[bs];
		Ciphertext* gsrot = new Ciphertext[gs];
		for(long i = 1; i < logrn; i++) {
			gap = n / (radix * pow(radix, i));
			tmpCipher = scheme.rightRotateFast(cipher, (radix - 1) * gap);
			bsrot[0] = tmpCipher;
			for(long j = 0; j < bs-1; j++) { // baby-step
				bsrot[j+1] = scheme.leftRotateFast(tmpCipher, (j+1) * gap);
			}
			for(long j = 0; j < gs; j++) {
				gsrot[j] = scheme.multByPoly(bsrot[0], encodeDFTpoly[i][j*bs], logp);
				for(long k = 1; k < bs; k++) {
					if((j * bs + k) < (2 * radix - 1)) {
						Ciphertext tmp = scheme.multByPoly(bsrot[k], encodeDFTpoly[i][j*bs+k], logp);
						scheme.addAndEqual(gsrot[j], tmp);
					}
				}
			}
			cipher = gsrot[0];
			for(long j = 0; j < gs-1; ++j) { // giant-step
				scheme.leftRotateFastAndEqual(gsrot[j+1], (j+1)*bs*gap);
			}
			for(long j = 1; j < gs; j++) {
				scheme.addAndEqual(cipher, gsrot[j]);
			}
			scheme.reScaleByAndEqual(cipher, logp);		
		}
		delete[] gsrot;
		delete[] bsrot;
	} else {
		for(long i = 1; i < logrn; i++) {
			gap = n / (radix * pow(radix, i));
			cipher = scheme.rightRotateFast(cipher, (radix - 1) * gap);
			tmpCipher = cipher;
			scheme.multByPolyAndEqual(cipher, encodeDFTpoly[i][0], logp);
			for(long j = 1; j < (2 * radix - 1); j++) {
				Ciphertext rot = scheme.leftRotateFast(tmpCipher, j * gap);
				scheme.multByPolyAndEqual(rot, encodeDFTpoly[i][j], logp);
				scheme.addAndEqual(cipher, rot);
			}
			scheme.reScaleByAndEqual(cipher, logp);
		}
	}
	return;
}

void DFTHelper::IDFT_RN(Ciphertext& cipher)
{
	// temperary things //
	ZZ* tmpPoly = new ZZ[ring.N]();
	Ciphertext tmpCipher;
	long bs, gs;

	// other step for invDFT //
	if(radix > 2) {
		// In this case, we use baby-step giant-step method //
		bs = 1 << (long)floor((log2(radix) + 1) / 2.);
		gs = (2 * radix) / bs;
		Ciphertext* bsrot = new Ciphertext[bs];
		Ciphertext* gsrot = new Ciphertext[gs];
		for(long i = 1; i < logrn; i++) {
			long gap = n / (radix * pow(radix, logrn-i));
			tmpCipher = scheme.rightRotateFast(cipher, (radix - 1) * gap);
			bsrot[0] = tmpCipher;
			for(long j = 1; j < bs; j++) { // baby-step
				bsrot[j] = scheme.leftRotateFast(tmpCipher, j * gap);
			}
			for(long j = 0; j < gs; j++) {
				gsrot[j] = scheme.multByPoly(bsrot[0], encodeinvDFTpoly[logrn-i][j*bs], logp);
				for(long k = 1; k < bs; k++) {
					if((j*bs+k) < (2*radix-1)) {
						Ciphertext tmp = scheme.multByPoly(bsrot[k], encodeinvDFTpoly[logrn-i][j*bs+k], logp);
						scheme.addAndEqual(gsrot[j], tmp);
					}
				}
			}
			cipher = gsrot[0];
			for(long j = 1; j < gs; j++) { // giant-step
				Ciphertext rot = scheme.leftRotateFast(gsrot[j], j*bs*gap);
				scheme.addAndEqual(cipher, rot);
			}
			scheme.reScaleByAndEqual(cipher, logp);		
		}
		delete[] gsrot;
		delete[] bsrot;
	} else {
		for(long i = 1; i < logrn; i++) {
			long gap = n / (radix * pow(radix, logrn-i));
			cipher = scheme.rightRotateFast(cipher, (radix - 1) * gap);
			tmpCipher = cipher;
			scheme.multByPolyAndEqual(cipher, encodeinvDFTpoly[logrn-i][0], logp);
			for(long j = 1; j < (2 * radix - 1); j++) {
				Ciphertext rot = scheme.leftRotateFast(tmpCipher, j * gap);
				scheme.multByPolyAndEqual(rot, encodeinvDFTpoly[logrn-i][j], logp);
				scheme.addAndEqual(cipher, rot);
			}
			scheme.reScaleByAndEqual(cipher, logp);
		}
	}

	// last step of invDFT //
	if(radix > 4) {
		// In this case, we use basy-step giant-step method
		// baby-step & giant-step method //
		bs = 1 << (long)floor((log2(radix)) / 2.);
		gs = (radix) / bs;
		Ciphertext* bsrotFirst = new Ciphertext[bs];
		Ciphertext* gsrotFirst = new Ciphertext[gs];		
		bsrotFirst[0] = cipher;
		for(long i = 1; i < bs; ++i) { // baby-step
			long gap = n / radix;
			bsrotFirst[i] = scheme.leftRotateFast(cipher, i * gap);
		}
		for(long i = 0; i < gs; ++i) {
			long gap = n / radix;
			gsrotFirst[i] = scheme.multByPoly(bsrotFirst[0], encodeinvDFTpoly[0][i*bs], logp);
			for(long j = 1; j < bs; ++j) {
				Ciphertext tmp = scheme.multByPoly(bsrotFirst[j], encodeinvDFTpoly[0][i*bs+j], logp);
				scheme.addAndEqual(gsrotFirst[i], tmp);
			}
		}
		cipher = gsrotFirst[0];
		for(long i = 1; i < gs; ++i) { // giant-step
			long gap = n / radix;
			Ciphertext rot = scheme.leftRotateFast(gsrotFirst[i], i*bs*gap);
			scheme.addAndEqual(cipher, rot);
		}
		delete[] bsrotFirst;
		delete[] gsrotFirst;
	} else {
		tmpCipher = cipher;
		scheme.multByPolyAndEqual(cipher, encodeinvDFTpoly[0][0], logp);
		for(long i = 1; i < radix; ++i) {
			long gap = n / radix;
			Ciphertext rot = scheme.leftRotateFast(tmpCipher, i * gap);
			scheme.multByPolyAndEqual(rot, encodeinvDFTpoly[0][i], logp);
			scheme.addAndEqual(cipher, rot);
		}
	}
	scheme.reScaleByAndEqual(cipher, logp);
	return;
}

DFTHelper::~DFTHelper() {
	delete[] DFTDecomp;
	delete[] invDFTDecomp;
	for(long i = 0; i < logrn; i++) {
		if(i == 0) {
			for(long j = 0; j < radix; j++) {
				delete[] encodeDFTpoly[i][j];
				delete[] encodeinvDFTpoly[i][j];
			}
		}
		else {
			for(long j = 0; j < 2 * radix - 1; j++) {
				delete[] encodeDFTpoly[i][j];
				delete[] encodeinvDFTpoly[i][j];
			}
		}
		delete[] encodeDFTpoly[i];
		delete[] encodeinvDFTpoly[i];
	}
	delete[] encodeDFTpoly;
	delete[] encodeinvDFTpoly;
}
