#include "stdafx.h"
#include "dftDll.h"
#include <stdexcept>

using namespace std;

namespace dft
{
	/**
	* Calculated Discrete Fourier Transform
	*  x - input array
	*  y - DFT result
	*/
	void Mydft::dft_forward(vector<double> &x, vector<complex<double>> &y)
	{
		int N = size(x);
		double t = -2 * acos((double)-1) / N;
		for (int k = 0;k < N;++k) {
			complex<double> sum(0.0, 0.0);
			for (int n = 0; n < N; ++n) {
				complex<double> expon(0.0, t*k*n);
				sum += (x[n]) * exp(expon);
			}
			y.push_back(sum);
		}
	}

	/**
	* Calculated Inverse Discrete Fourier Transform
	*  y - input array
	*  x - IDFT result
	*/
	void Mydft::dft_inverse(vector<complex<double>> &y, vector<complex<double>> &x)
	{
		double N = size(y);
		double t = 2 * acos((double)-1) / N;
		for (int n = 0;n < N;++n) {
			complex<double> sum(0.0, 0.0);
			for (int k = 0; k < N; ++k) {
				complex<double> expon(0.0, t*k*n);
				sum += y[k] * exp(expon);
			}
			x.push_back(sum / complex<double>(N, 0));
		}
	}

	/**
	* Calculated Analytic Signal of input signal x using DFT
	*  x - input array
	*  y - analytic signal
	*/
	void Mydft::analytic_signal_dft(vector<double> &x, vector<complex<double>> &y) {
		int N = size(x);
		vector<complex<double>> x_dft;
		dft_forward(x, x_dft);
		//	x_dft= DFT of x
		for (int i = 1;i < N / 2;++i) {
			x_dft[i] = complex<double>(2, 0)*x_dft[i]; x_dft[i + N / 2] = 0;
		}
		//  now x_dft contains only the non-negative frequency components of x_dft
		dft_inverse(x_dft, y);
		//  y= Analytic signal of input signal x
	}
}

