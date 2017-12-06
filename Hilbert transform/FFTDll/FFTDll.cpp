// FFTDll.cpp: определяет экспортированные функции для приложения DLL.
//

#include "stdafx.h"
#include "fftDllh.h"
#include <stdexcept>

using namespace std;

namespace fft
{
	/**
	* Calculated Fast Fourier Transform if N=2^k
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

