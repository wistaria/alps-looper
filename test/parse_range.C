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

#include <alps/expression.h>
#include <alps/parameters.h>
#include <boost/spirit/core.hpp>
// #include <boost/spirit/actor.hpp>
// #include <boost/spirit/utility/lists.hpp>
#include <cstdlib>
#include <iostream>
#include <limits>

template<class T>
bool parse_range(char const* str, std::pair<T, T>& r,
                 const alps::Parameters& params = alps::Parameters(),
                 T min_default = std::numeric_limits<T>::min(),
                 T max_default = std::numeric_limits<T>::max())
{
  typedef T value_type;
  using namespace boost::spirit;

  std::string min_str, max_str;
  bool success = parse(str,
    '[' >> (*(anychar_p - ':' - ']'))[assign_a(min_str)] >> ':' >> 
    (*(anychar_p - ':' - ']'))[assign_a(max_str)] >> ']',
    space_p).full;
  if (!success) return false;

  if (min_str.empty()) {
    r.first = min_default;
  } else {
    if (!alps::can_evaluate(min_str, params)) return false;
    r.first = alps::evaluate<double>(min_str, params);
  }

  if (max_str.empty()) {
    r.first = max_default;
  } else {
    if (!alps::can_evaluate(max_str, params)) return false;
    r.first = alps::evaluate<double>(max_str, params);
  }

  return true;
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::string str;
  while (std::getline(std::cin, str)) {
    if (str[0] == 'q') break;

    std::pair<double, double> r;
    if (parse_range(str.c_str(), r)) {
      std::cout << "parsing succeeded: " << str << std::endl;
      std::cout << "\tmin : " << r.first << std::endl
                << "\tmax : " << r.second << std::endl;
    } else {
      std::cout << "parsing failed: " << str << std::endl;
    }
  }

  while (std::getline(std::cin, str)) {
    if (str[0] == 'q') break;

    std::pair<unsigned int, unsigned int> r;
    if (parse_range(str.c_str(), r)) {
      std::cout << "parsing succeeded: " << str << std::endl;
      std::cout << "\tmin : " << r.first << std::endl
                << "\tmax : " << r.second << std::endl;
    } else {
      std::cout << "parsing failed: " << str << std::endl;
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
}
