#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  if (argc == 3) {
    int L = boost::lexical_cast<std::size_t>(argv[1]);
    int count = boost::lexical_cast<std::size_t>(argv[2]);
    std::vector<double> vec;
    vec.reserve(L);
    // std::vector<double> vec(L, 0);
    boost::timer timer;
    for (int c = 0; c < count; ++c ) {
      vec.resize(0);
      for (int i = 0; i < L; ++i)
        vec.push_back(0.1);
    }
    double elp = timer.elapsed();
    std::cout << "vec size = " << vec.size() << std::endl
              << "count = " << count << std::endl
              << "elapsed = " << elp << " sec\n"
              << "speed = " << 1.0 * L * count / elp << " elements/sec\n";
  } else {
    std::cerr << "usage: " << argv[0] << " L count\n";
    return 127;
  }
}
