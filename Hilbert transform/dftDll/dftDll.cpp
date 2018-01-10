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
	
	
/**
 * Filter input signal
 * x - input array
 * X - filter result
 * w0 - lower cutoff frequency
 * w1 - higher cutoff frequency
 * Ep - passband feature for lpf(smaller is better)
 * Es - stopband feature for lpf(larger is better)
 * fil = 1 : Butterworth approximation; fil = 2 : type first Chebyshev approximation
 * pass = 1 : low-pass filter; pass = 2 : high-pass filter
 */
void Mydft::some_filter(vector<double> x, vector<complex<double>> &X, double w0, double w1, double Ep, double Es, int fil, int pass) {
    if(w1 > x.size()/2) {
        cout<<" Sampling rate is not enought ";
        return;
    }
    vector<complex<double>> z,H;
    switch(fil){
        case 1:{ //count filter with Butterworth approximation
            int N = log(Es / Ep) / log(w1 / w0) + 1; //count number of poles in the filter
            for (int i = 0; i < int(x.size())/2 ; ++i) {
                H.push_back((1 / sqrt(1 + pow(((i * 2 * M_PI) / w0), 2 * N) * pow(Ep, 2))));
            }
            break;
        }
        case 2:{ // count filter with type first Chebyshev approximation
            int N = log((Es / Ep)-sqrt(pow((Es/Ep), 2) - 1)) / log((w1 / w0)-sqrt(pow((w1/w0), 2) - 1)) + 1; //count number of poles in the filter
            for (int i = 0; i < int(x.size())/2 ; ++i) {
                complex<double> m;
                if(i * 2 * M_PI > w0){
                    m = complex<double>(-1, 0) * log(complex<double>(double(i * 2 *M_PI) / w0, 0) - complex<double>( sqrt(-1 + pow(double(i * 2 *M_PI) / w0, 2)), 0)) * complex<double>(0, 1);
                }
                else {
                    m = log(complex<double>(double(i * 2 *M_PI) / w0, 0) + complex<double>( sqrt(1 - pow(double(i * 2 *M_PI) / w0, 2)), 0) * complex<double>(0, 1)) * complex<double>(0, -1);
                }
                
                H.push_back((complex<double>(1, 0) / sqrt(complex<double>(1, 0) + pow(cos(complex<double>(N, 0) * m), complex<double>(2, 0)) * pow(Ep, 2))));
            }
            break;
        }
    }
    
    switch(pass){
        case 1:{ // lpf
            break;
        }
        case 2:{ // hpf
            for (int i = 0; i < int(x.size())/2 ; ++i) {
                H[i] = complex<double>(1, 0)-H[i];
            }
            break;
        }
    }
    
    dft_forward(x, z);
    for (int i = 0; i < int(z.size())/2 ; ++i) {
        z[i] *= H[i];
    }
    for(int i = int(z.size())/2; i< z.size(); ++i) {
        z[i] = -z[z.size()-i];
    }
    dft_inverse(z, X);
}
	

}

