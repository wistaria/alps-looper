#include "vectorhelper.h"
#include <boost/random.hpp>
#include <iostream>

template<class C>
void check(const C& c)
{
  looper::index_helper<C> helper(c);
  for (int i = 0; i < c.size(); ++i) assert(i == helper.index(&c[i]));
}

int main()
{
  boost::mt19937 rng;
  boost::uniform_int<> uniform_int(1, 10);

  std::vector<int> a(100);
  check(a);
  std::cout << "std::vector(" << a.size() << ") OK\n";

  std::deque<int> b(uniform_int(rng));
  check(b);
  for (;;) {
    int n = uniform_int(rng);
    if (rng() < 0.5) {
      for (int i = 0; i < n; ++i) b.push_front(0);
    } else {
      for (int i = 0; i < n; ++i) b.push_back(0);
    }
    check(b);
    std::cout << "std::deque(" << b.size() << ") OK\n";
    if (b.size() > 300) break;
  }
  assert(looper::index_helper<>::index(b, &b[b.size() / 2]) == b.size() / 2);
}
