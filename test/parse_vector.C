/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2005-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <boost/classic_spirit.hpp>
#include <iostream>
#include <cstdlib>

template<class T>
bool parse_vector(char const* str, T& v)
{
  using namespace boost::spirit;

  subrule<0> vec;
  subrule<1> elem;

  return parse(str,
    (
    vec  = '[' >> elem >> *(',' >> elem) >> ']',
    elem = (+(anychar_p - ',' - ']'))[push_back_a(v)]
    ),
    space_p).full;
}

template<class T>
bool parse_vectors(char const* str, T& c)
{
  using namespace boost::spirit;

  typename T::value_type v;
  v.clear();

  subrule<0> cont;
  subrule<1> vec;
  subrule<2> elem;

  return parse(str,
    (
    cont  = '[' >> vec >> *(',' >> vec) >> ']',
    vec   = eps_p[clear_a(v)] >>
            '[' >> elem >> *(',' >> elem) >> ']' >>
            eps_p[push_back_a(c, v)],
    elem  = (+(anychar_p - ',' - ']'))[push_back_a(v)]
    ),
    space_p).full;
}

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::string str;

  while (std::getline(std::cin, str)) {
    if (str[0] == 'q') break;

    std::vector<std::string> v;
    if (parse_vector(str.c_str(), v)) {
      std::cout << "parsing succeeded: " << str << std::endl;
      for (int i = 0; i < v.size(); ++i)
        std::cout << '\t' << i << ": " << v[i] << std::endl;
    } else {
      std::cout << "parsing failed: " << str << std::endl;
    }
  }

  while (std::getline(std::cin, str)) {
    if (str[0] == 'q') break;

    std::vector<std::vector<std::string> > c;
    if (parse_vectors(str.c_str(), c)) {
      std::cout << "parsing succeeded: " << str << std::endl;
      for (int i = 0; i < c.size(); ++i)
        for (int j = 0; j < c[i].size(); ++j)
          std::cout << '\t' << i << ' ' << j << ": " << c[i][j] << std::endl;
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
