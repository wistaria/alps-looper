#include <looper/atomic.h>
#include <looper/atomic_impl.h>
#include <iostream>
#include <vector>

int main() {
  int x = 0;
  std::cout << x << std::endl;
  bool res0 = looper::compare_and_swap(x, 0, 3);
  std::cout << x << ' ' << res0 << std::endl;
  bool res1 = looper::compare_and_swap(x, 0, 5);
  std::cout << x << ' ' << res1 << std::endl;
  bool res2 = looper::compare_and_swap(x, 3, 5);
  std::cout << x << ' ' << res2 << std::endl;

  const int n = 10000000;
  looper::atomic_counter counter;
  std::vector<int> check(n, 0);
  #pragma omp parallel for
  for (int i = 0; i < n; ++i) {
    check[counter++] = 1;
    // check[counter.non_atomic_increment()] = 1;
  }
  std::cout << counter() << std::endl;
  for (int i = 0; i < n; ++i) {
    if (check[i] != 1) {
      std::cout << "Error\n";
      break;
    }
  }
}
