/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004-2007 by Synge Todo <wistaria@comp-phys.org>
*
* This software is published under the ALPS Application License; you
* can use, redistribute it and/or modify it under the terms of the
* license, either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License
* along with this software; see the file LICENSE. If not, the license
* is also available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include <looper/poisson_distribution.h>
#include <alps/parameter.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random.hpp>
#include <vector>

int main() {
  alps::Parameters params(std::cin);

  double mean = static_cast<int>(params["MEAN"]);
  int count = params["COUNT"];

  //
  // random number generator
  //

  boost::mt19937 engine;
  boost::variate_generator<boost::mt19937, looper::poisson_distribution<> >
    rng(engine, looper::poisson_distribution<>(mean));

  //
  // test
  //

  std::vector<int> bins(int(5 * mean), 0);
  boost::posix_time::ptime start = boost::posix_time::microsec_clock::universal_time();
  for (int c = 0; c < count; ++c) {
    int r = rng();
    if (r < 5 * mean) ++bins[r];
  }
  boost::posix_time::ptime stop = boost::posix_time::microsec_clock::universal_time();
  double time = 1.e-6 * (stop - start).total_microseconds();
  std::cerr << "count = " << count
            << ", total = " << time << " sec, "
            << time / count << " sec/sweep\n";

  std::cout << std::setprecision(3);
  double poi = std::exp(-mean);
  for (int r = 0; r < bins.size(); ++r) {
    std::cout << r << ' ' << poi << ' '
              << (double)bins[r] / count << ' '
              << std::sqrt((double)bins[r]) / count << std::endl;
    poi *= mean / (r+1);
  }
}
