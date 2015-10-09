#include <boost/timer.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <vector>
#include <typeinfo>

// #include <fj_tool/fjsamp.h>

struct pod_t {
  pod_t& operator+=(pod_t const& p) {
    a += p.a;
    b += p.b;
    c += p.c;
    d += p.d;
    return *this;
  }
  int a;
  float b;
  double c;
  unsigned int d;
};

struct non_pod_t {
  non_pod_t() : a(0), b(0), c(0), d(0) {}
  non_pod_t& operator+=(non_pod_t const& np) {
    a += np.a;
    b += np.b;
    c += np.c;
    d += np.d;
    return *this;
  }
  int a;
  float b;
  double c;
  unsigned int d;
};

template<class POD>
void check_pod (POD p_ob) {
  if (p_ob.a != 0 || p_ob.b != 0 || p_ob.c != 0 || p_ob.d) {
    // POD object is NOT initialized by constructor
    std::cout << typeid(POD).name() << " is POD type" << std::endl;
  } else {
    // Non POD object is initialized by constructor
    std::cout << typeid(POD).name() << " is NOT POD type!!" << std::endl;
  }
}


void output(const int size, const int count, const int L, double elp) {
  std::cout << "vec size = " << size << std::endl
            << "count = " << count << std::endl
            << "elapsed = " << elp << std::endl
            << "speed = " << 1.0 * L * count / elp << " elements/sec" << std::endl;
}

int main(int argc, char** argv)
{
  // POD object doesn't call constructor 
  pod_t ob;
  non_pod_t nob;
  check_pod(ob);
  check_pod(nob);

  // Initialize POD object manually
  ob = pod_t();
  
  if (argc == 3) {
    int L = boost::lexical_cast<std::size_t>(argv[1]);
    int count = boost::lexical_cast<std::size_t>(argv[2]);
    std::vector<pod_t> vec_pod(L);
    std::vector<non_pod_t> vec_npod(L);
    boost::timer timer;
    // fpcoll_start();
    for (int c = 0; c < count; ++c) {
      // vec_pod.resize(L);
      for (int i = 0; i < L; ++i)
        vec_pod[i] = ob;
    }
    // fpcoll_stop();
    double elp_pod = timer.elapsed();

    timer.restart();
    // fpcoll_start();
    for (int c = 0; c < count; ++c) {
      // vec_npod.resize(L);
      for (int i = 0; i < L; ++i)
        vec_npod[i] = nob;
    }
    // fpcoll_stop();
    double elp_npod = timer.elapsed();

    std::cout << "pod:" << std::endl;
    output(vec_pod.size(), count, L, elp_pod);
    std::cout << "non_pod:" << std::endl;
    output(vec_npod.size(), count, L, elp_npod);
  } else {
    std::cerr << "usage: " << argv[0] << " L count" << std::endl;
    std::exit(127);
  }
  return 0;  
}
