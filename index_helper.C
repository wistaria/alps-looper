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
  boost::uniform_int<> uniform_int(1, 5);

  std::vector<int> a(uniform_int(rng));
  check(a);
  for (;;) {
    a.push_back(uniform_int(rng));
    check(a);
    if (a.size() > 1000) break;
  }
  assert(looper::index_helper<>::index(a, &a[a.size() / 2]) == a.size() / 2);

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
    if (b.size() > 1000) break;
  }
  assert(looper::index_helper<>::index(b, &b[b.size() / 2]) == b.size() / 2);

  std::cout << "OK\n";
}
