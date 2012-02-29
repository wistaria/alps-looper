#include <iostream>
#include <vector>
#include <boost/lexical_cast.hpp>

#ifdef _OPENMP
# include <omp.h>
#else
inline int omp_get_max_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
#endif

int main(int argc, char** argv) {
  int n = (argc >= 2) ? boost::lexical_cast<std::size_t>(argv[1]) : 10;
  std::vector<int> v(n);
  std::cerr << "number of threads = " << omp_get_max_threads() << std::endl;
  #pragma omp parallel
  {
    int tid = omp_get_thread_num();
    #pragma omp for schedule(static)
    for (int i = 0; i < n; ++i) v[i] = tid;
  }
  for (int i = 0; i < n; ++i) std::cerr << i << ' ' << v[i] << std::endl;
}
