#ifndef ALPS_ENABLE_TIMER
# define ALPS_ENABLE_TIMER
#endif

#include "looper/timer.hpp"
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <vector>
#ifdef _OPENMP
# include <omp.h>
#endif


int main(int argc, char** argv) {
  int counts = (argc >= 2) ? boost::lexical_cast<std::size_t>(argv[1]) : 1 << 10;
  int pad_min = 0;
  int pad_max = 16;
  
#ifdef _OPENMP
  const int num_threads = omp_get_max_threads();
  const int thread_id = omp_get_thread_num();
#else
  const int num_threads = 1;
  const int thread_id = 0;
#endif
  std::cerr << "number of threads = " << num_threads << std::endl;
  std::cerr << "number of counts = " << counts << std::endl;

  alps::parapack::timer timer;
  timer.registrate(0, "timer");
  
  // for (int p = pad_min; p <= pad_max; ++p) {
  for (int q = pad_min; q <= pad_max; ++q) {
    int p = 1 << q;
    std::vector<double> vec(num_threads * (p + 1));
    timer.clear();
    timer.start(0);
    #pragma omp parallel
    {
      #ifdef _OPENMP
      int index = omp_get_thread_num() * (p + 1);
      #else
      int index = 0;
      #endif
      for (int i = 0; i < counts; ++i) {
        vec[index] += i;
      }
    }
    timer.stop(0);
    std::cerr << "padding = " << p * sizeof(double) << " byte: "
              << "elapse = " << timer.get_measurements()[0] << " sec\n";
    // for (int i = 0; i < vec.size(); ++i) std::cerr << i << ' ' << vec[i] << std::endl;
  }
}
