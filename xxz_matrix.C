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
	      << std::endl
	      << xxz.matrix() << std::endl;

    boost::tuple<double, double, double> p =
      looper::fit_xxz_matrix(s0, s1, xxz.matrix());
    std::cout << "fitting result: e0 = " << p.get<0>()
	      << ", Jxy = " << p.get<1>()
	      << ", Jz = " << p.get<2>()
	      << std::endl;
  }
}
