#include "xxz_matrix.h"
#include <iostream>

int main()
{
  int n;
  std::cin >> n;
  for (int i = 0; i < n; ++i) {
    double s0_in, s1_in, e0, jxy, jz;
    std::cin >> s0_in >> s1_in >> e0 >> jxy >> jz;
    alps::half_integer<int> s0(s0_in);
    alps::half_integer<int> s1(s1_in);
    looper::xxz_matrix<> xxz(s0, s1, e0, jxy, jz);

    std::cout << "input parameters: S0 = " << s0 << ", S1 = " << s1
	      << ", e0 = " << e0 << ", Jxy = " << jxy << ", Jz = " << jz
	      << std::endl << xxz << std::endl;

    boost::tuple<bool, double, double, double> fit =
      looper::fit_xxz_matrix(s0, s1, xxz);

    if (fit.get<0>()) {
      std::cout << "fitting result: e0 = " << fit.get<1>()
		<< ", Jxy = " << fit.get<2>()
		<< ", Jz = " << fit.get<3>()
		<< std::endl;
    } else {
      std::cout << "fitting failed\n";
    }
  }
}
