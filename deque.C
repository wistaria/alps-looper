#include "deque.h"
#include <iostream>

int main()
{
  std::deque<int> a(10);
  a.push_front(3);

  looper::deque_index_helper<std::deque<int> > h(a);
  for (int i = 0; i < a.size(); ++i) {
    std::cout << i << ' ' << h.index(&a[i]) << std::endl;
  }
}
