#include "deque.h"
#include <iostream>

template<class C>
void check(const C& c)
{
  looper::index_helper<C> helper(c);
  bool check = true;
  for (int i = 0; i < c.size(); ++i) {
    if (i != helper.index(&c[i])) {
      check = false;
      break;
    }
  }

  if (!check) {
    std::cout << "ERROR occured\n";
    helper.output();
    for (int i = 0; i < c.size(); ++i) 
      std::cout << i << ' ' << helper.index(&c[i], true) << ' ' << &c[i]
		<< std::endl;
    std::abort();
  }

  std::cout << c.size() << ' ' << helper.num_chunks() << std::endl;
}

int main()
{
  std::deque<int> a(30);
  check(a);
  for (;;) {
    a.push_front(3);
    a.push_front(3);
    a.push_front(3);
    check(a);
    a.push_back(3);
    a.push_back(3);
    check(a);
    if (a.size() > 300) break;
  }
}
