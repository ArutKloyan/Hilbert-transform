#ifdef MATHFUNCSDLL_EXPORTS
#define MATHFUNCSDLL_API __declspec(dllexport) 
#else
#define MATHFNCSDLL_API __declspec(dllimport) 
#endif
#include <complex>
#include <vector>
using namespace std;
namespace dft
{
	class Mydft
	{
	public:
		static MATHFNCSDLL_API void dft_forward(vector<double> &x, vector<complex<double>> &y);
		static MATHFNCSDLL_API void dft_inverse(vector<complex<double>> &y, vector<complex<double>> &x);
		static MATHFNCSDLL_API void analytic_signal_dft(vector<double> &x, vector<complex<double>> &y);
		static MATHFNCSDLL_API void some_filter(vector<double> x, vector<complex<double>> &X, double w0, double w1, double Ep, double Es, int fil, int pass);
	};
}
