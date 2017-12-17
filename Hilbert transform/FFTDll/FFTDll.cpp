#include "stdafx.h"
#include "fftDllh.h"
#include <stdexcept>

using namespace std;

namespace fft
{
	/**
	* Calculated Fast Fourier Transform if N=2^k
	*  y - input array and FFT result
	*  N - number of points
	*/
	void Myfft::fft_forward(vector<complex<double>> &y) {
		int N=size(y);
	vector<complex<double>> b;
	for (int i = 0; i < N; ++i) {
		b.push_back(y[i]);
	}
	for (int j = N; j != 2; j /= 2) {
		for (int k = 0; k*j < N; ++k) {
			for (int i = 0; i < j / 2; ++i) {
				b[k*j + i] = y[k*j + 2 * i];
				b[(k*j + i) + j / 2] = y[k*j + 2 * i + 1];
			}
			for (int i = 0; i < j; ++i) {
				y[(k*j + i)] = b[(k*j + i)];
			}
		}
	}
	
	complex<double> w1, w2;
	for (int j = 2; j <= N; j *= 2) {
		for (int k = 0; k*j < N; ++k) {
			for (int i = 0; i < j / 2; ++i) {
				complex<double> expon = exp(complex<double>(0, -2.*(3.1415926535)*i / j));
				w1 = y[k*j + i];
				w2 = y[k*j + i + j / 2] * expon;
				y[k*j + i] = w1 + w2;
				y[k*j + i + j / 2] = w1 - w2;
			}
		}
	}
}
	
	
	
	
	
	
	/**
	* Calculated Fast Fourier Transform with recursion if N=2^k
	*  X - input array and FFT result
	*  N - number of points
	*/
	void Myfft::fft_forward_ct(complex<double> * X, int N) {
		if (N < 2) {
			return;
		}
		vector<complex<double>> b;
		for (int i = 0; i < N / 2; i++) {
			b.push_back(X[2 * i + 1]);
		}
		for (int i = 0; i < N / 2; i++) {
			X[i] = X[2 * i];
		}
		for (int i = 0; i < N / 2; i++) {
			X[i + N / 2] = b[i];
		}
		fft_forward_ct(X, N / 2);   // recurse even items
		fft_forward_ct(&(X[N / 2]), N / 2);   // recurse odd  items
		// combine results of two recursions
		for (int k = 0; k < N / 2; k++) {
			complex<double> e = X[k];   // even
			complex<double> o = X[k + N / 2];   // odd											
			complex<double> expon = exp(complex<double>(0, -2.*(3.1415926535)*k / N));
			X[k] = e + expon * o;
			X[k + N / 2] = e - expon * o;
		}
	}

	/**
	* Calculated Inverse Fast Fourier Transform
	*  y - input array and IDFT result
	*/
	void Myfft::fft_inverse_ct(vector<complex<double>>& y) {
		int N = size(y);
		for (int i = 0; i < N; ++i) {
			y[i] = complex<double>(imag(y[i]), real(y[i])); 
		}
		fft_forward_ct(&(y[0]), N);
		for (int i = 0; i < N; ++i) { 
			y[i] = complex<double>(imag(y[i]) / N, real(y[i]) / N); 
		}
	}

	/**
	* Calculated Analytic Signal of input signal x using FFT
	*  x - input array
	*  y - analytic signal
	*/
	void Myfft::analytic_signal_fft(vector<double> &x, vector<complex<double>> &y) {
		int N = size(x);
		for (int i = 0; i < size(x); ++i) y.push_back(complex<double>(x[i], 0.0));
		fft_forward_ct(&(y[0]), N);
		//	y= FFT of x
		for (int i = 1;i < N / 2;++i) {
			y[i] = complex<double>(2, 0)*y[i]; y[i + N / 2] = 0;
		}
		//  now y contains only the non-negative frequency components of x_dft
		fft_inverse_ct(y);
		//  y= Analytic signal of input signal x
	}
	
	
	/**
	* Calculated multiple of two polynomials using FFT
	*  x - first polynomial
	*  y - second polynomial
	*  mul - multiple of x and y
	*/
	void Myfft::polinomials_mul_fft(vector<double> x, vector<double> y, vector<double> & mul) {
	vector<complex<double>> x1, y1, mul1;
	int N1, N2;
	N1 = size(x);
	N2 = size(y);
	for (int i = 0; i < N1; ++i) {
		x1.push_back(complex<double>(x[i], 0.0));
	}
	for (int i = 0; i < N2; ++i) {
		y1.push_back(complex<double>(y[i],0.0));
	}
		
	//degree of multiple is a sum of degrees
	int N = N1 + N2;
		
	//FFT needs a degree of 2
	//find the power of two which is closest to N
	int i = 2;
	while (i < N) {
		i <<= 1;
	}
		
	N = i;
	for (int i = N1; i < N; ++i) {
		x1.push_back(complex<double>(0, 0));
	}
	for (int i = N2; i < N; ++i) {
		y1.push_back(complex<double>(0, 0));
	}
	fft_forward_ct(&(x1[0]), N);
	fft_forward_ct(&(y1[0]), N);
	
	for (int i = 0; i < N; ++i) {
		mul1.push_back(x1[i] * y1[i]);
	}
		fft_inverse_ct(mul1);
	for (int i = 0; i < N; ++i)  mul.push_back(real(mul1[i]));
	
}


	/**
	* if the difference between x and y are larger than test,then programm is working bad
	* if the difference between x and y are less than test,then programm is working good
	*  x - dft result
	*  y - fft result
	*  test - maximum difference between elements of x and y
	*/
	void Myfft::fft_dft(vector<double> &x, vector<complex<double>> &y, double test) {
		int N = size(x);
		for (int i = 0; i < N; ++i) {
			if (((x[i] - real(y[i])) > test) || ((real(y[i]) - x[i]) > test)) { 
//				std::cout << "Something went wrong" << endl; return;
			}
		}
//				std::cout << "Everything is OK" << endl;
	}
}

