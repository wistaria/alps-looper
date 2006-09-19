/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2006 by Synge Todo <wistaria@comp-phys.org>
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

#include <boost/regex.hpp>
#include <iostream>
#include <stdexcept>

void match(std::string const& str, std::string const& ex)
{
  std::cout << "str = \"" << str << "\", ex = \"" << ex << "\", ";
  if (regex_match(str, boost::regex(ex)))
    std::cout << "result = matched\n";
  else
    std::cout << "result = not matched\n";
}

void replace(std::string const& str, std::string const& ex,
             std::string const& fmt)
{
  std::cout << "str = \"" << str << "\", ex = \"" << ex << "\", fmt = \""
            << fmt << "\", ";
  std::string out =
    regex_replace(str, boost::regex(ex), fmt,
                  boost::match_default | boost::format_all);
  std::cout << "result = " << out << std::endl;
}

int main(int , char**)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  // match test
  match("SSE QWL", "(.*)");
  match("SSE QWL", "QWL");
  match("SSE QWL", "(.*)QWL$");
  match("SSE QWL", "^(.*)QWL$");
  match("SSE QWL", "(.*)QWL(.*)");
  match("SSE QWL A", "(.*)QWL$");
  match("SSE QWL A", "(.*)QWL(.*)$");
  match("SSE QWL A", "QWL");

  // replace test
  replace("a.out.xml", "out.xml$", "ox.bak");
  replace("a.out.xmla", "out.xml$", "ox.bak");
  replace("a.out.xmla", "out.xml", "ox.bak");
  replace("a.outxxml", "out.xml", "ox.bak");
  replace("a.outxxml", "out\\.xml", "ox.bak");
  replace("a.out.xml", "out\\.xml", "ox.bak");
  replace("a.outxxml", "out[:punkt:]xml", "ox.bak");
  replace("a.out.xml", "out[:punkt:]xml", "ox.bak");
  replace("a.out.xml", "\\.out\\.xml$", "");
  replace("Energy Density", "[A-Z]", "\\l$&");
  replace("Free Energy Density", "\\s", "_");
  replace("Free Energy Density", "\\sDensity$", "");
  replace("Free Energy Density A", "\\sDensity$", "");
  replace("Free Energy Density", "(\\sDensity$)|([A-Z])|(\\s)",
          "(?2\\l$2)(?3_)");
  replace("Staggered Susceptibility", "(\\sDensity$)|([A-Z])|(\\s)",
          "(?2\\l$2)(?3_)");

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}
