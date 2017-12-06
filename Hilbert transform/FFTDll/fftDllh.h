#ifdef MATHFUNCSDLL_EXPORTS
#define MATHFUNCSDLL_API __declspec(dllexport) 
#else
#define MATHFUNCSDLL_API __declspec(dllimport) 
#endif
#include <complex>
#include <vector>
using namespace std;
namespace fft
{
	class Myfft
	{
	public:
		static MATHFUNCSDLL_API void fft_forward_ct(complex<double> * X, int N);
		static MATHFUNCSDLL_API void fft_inverse_ct(vector<complex<double>>& y);
		static MATHFUNCSDLL_API void analytic_signal_fft(vector<double> &x, vector<complex<double>> &y);
		static MATHFUNCSDLL_API void fft_dft(vector<double> &x, vector<complex<double>> &y, double test);
	};
}