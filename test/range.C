/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2005 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/util.h>
// #include <alps/expression.h>
// #include <alps/parameters.h>
// #include <boost/spirit/core.hpp>
// #include <boost/spirit/actor.hpp>
// #include <boost/spirit/utility/lists.hpp>
#include <cstdlib>
#include <iostream>
#include <limits>

int main()
{
  std::string str;
  while (std::getline(std::cin, str)) {
    if (str[0] == 'q') break;
    std::cout << "parse " << str << ": ";
    try {
      looper::integer_range<int> r(str);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }
  while (std::getline(std::cin, str)) {
    if (str[0] == 'q') break;
    std::cout << "parse " << str << ": ";
    try {
      looper::integer_range<unsigned int> r(str);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }
}
