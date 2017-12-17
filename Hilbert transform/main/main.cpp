#include "stdafx.h"
#include "fftDllh.h"
#include "dftDll.h"

double input_func(double);
void generate_test_input_signal(vector<double>&);
void generate_polynomial(vector<double>&);
void check_result(vector<double>&, vector<complex<double>>&, double);

int main()
{
	vector<double> x;
	vector<complex<double>> y, z;
	generate_test_input_signal(x);
	//  x = input signal
	fft::Myfft::analytic_signal_fft(x, y);
	//  y = analytic signal of x using fft
	dft::Mydft::analytic_signal_dft(x, z);
	//  y = analytic signal of x using dft
	check_result(x, y, 0.00001);
	check_result(x, z, 0.00001);
	system("pause");
	return 0;
	return 0;
}


/**
* Calculated element of input signal
*  x - input unit
*/
double input_func(double x) {
	return 5 * sin(3 * x) + 13 * sin(4 * x) + 4 * sin(11 * x);
}


/**
* Calculated input signal
*  x - input array
*/
void generate_test_input_signal(vector<double> &x) {
	int N=16;
//	printf("The number of values = ");
//	cin >> N;
	printf("\n");
	double t = 2 * acos((double)-1) / N;
	for (int i = 0;i < N; ++i) {
		x.push_back((input_func(t*i)));
	}
}


/**
* Calculated input polynomial
*  x - input array
*/
void polynomial(vector<double> & x)
{
	int N;
	cin >> N;
	for (int i = 0; i <= N; ++i) {
		double a;
		cin >> a;
		x.push_back(a);
	}
}

/**
* if the difference between x and y are larger than epsilon,then programm is working bad
* if the difference between x and y are less than epsilon,then programm is working good
* the real and imaginary part of the analytic signal should satisfy the property of orthogonality
*  x - input array
*  y - analytic signal
*  test - maximum difference between elements of x and y
*/
void check_result(vector<double> &x, vector<complex<double>> &y, double epsilon) {
	int N = size(x);
	for (int i = 0; i < N; ++i) {
		if (((x[i] - real(y[i])) > epsilon) || ((real(y[i]) - x[i]) > epsilon)) {
			printf("The difference between x and y are larger than epsilon\n"); return;
		}
	}
	complex<double> sum(0.0, 0.0);
	for (int i = 0; i < N; ++i) {
		sum += (real(y[i])*imag(y[i]));
	}
	if (((real(sum)) > epsilon) || ((imag(sum)) > epsilon)) {
		printf("The real and imaginary part of the analytic signal do not satisfy the property of orthogonality\n"); return;
	}
	printf("Everything is OK\n");
}

