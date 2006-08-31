/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2005-2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <looper/integer_range.h>
#include <cstdlib>
#include <iostream>
#include <string>

int main()
{
  alps::Parameters params;
  std::cin >> params;
  std::cout << "Parameters:\n" << params;

  std::string str;
  std::cout << "Test for integer_range<int>:\n";
  while (std::getline(std::cin, str)) {
    if (str[0] == '\0') break;
    std::cout << "parse " << str << ": ";
    try {
      looper::integer_range<int> r(str, params);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }
  std::cout << "Test for integer_range<unsigned int>:\n";
  while (std::getline(std::cin, str)) {
    if (!std::cin || str[0] == '\0') break;
    std::cout << "parse " << str << ": ";
    try {
      looper::integer_range<unsigned int> r(str, params);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }

  looper::integer_range<int> r(0, 5);
  std::cout << "initial: " << r << std::endl;
  std::cout << "7 is included? " << r.is_included(7) << std::endl;
  r = 3;
  std::cout << "3 is assigned: " << r << std::endl;
  r.include(8);
  std::cout << "8 is included: " << r << std::endl;
  std::cout << "7 is included? " << r.is_included(7) << std::endl;
}
