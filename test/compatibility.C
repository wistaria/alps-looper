/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2005 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/node.h>
#include <iostream>

int main() {

#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << "[bond graph]\n";
  for (int g = 1; g < 5; ++g)
    for (int c0 = 0; c0 < 2; ++c0)
      for (int c1 = 0; c1 < 2; ++c1)
        std::cout << "Delta(g=" << g
                  << ", c0=" << (c0 == 0 ? '+' : '-')
                  << ", c1=" << (c1 == 0 ? '+' : '-')
                  << ") = "
                  << (looper::is_compatible(looper::bond_graph(0, g), c0, c1) ?
                      1 : 0)
                  << std::endl;

  std::cout << "[site graph]\n";
  for (int g = 1; g < 4; ++g)
    for (int c = 0; c < 2; ++c)
      std::cout << "Delta(g=" << g
                << ", c=" << (c == 0 ? '+' : '-')
                << ") = "
                << (looper::is_compatible(looper::site_graph(0, g), c) ? 1 : 0)
                << std::endl;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (const std::exception& excp) {
  std::cerr << excp.what() << std::endl;
  std::exit(-1); }
catch (...) {
  std::cerr << "Unknown exception occurred!" << std::endl;
  std::exit(-1); }
#endif

}