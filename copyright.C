/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: copyright.C 604 2004-01-16 08:35:21Z wistaria $
// copyright - print copyright and/or license information

#include "looper/copyright.h"
#include <boost/throw_exception.hpp>
#include <iostream>
#include <stdexcept>

struct options
{
  bool license;

  options(int argc, char *argv[]) : license(false) { parse(argc, argv); }

  void usage(int status, std::ostream& os = std::cerr) const
  {
    os << "NAME\n"
       << "     copyright - print copyright and/or license information\n\n"
       << "SYNOPSIS\n"
       << "     copyright [-l] [-h]\n\n"
       << "DESCRIPTION\n"
       << "     The options are as follows:\n\n"
       << "     -l          print license details.\n"
       << "     -h          print this help message and exit.\n\n";
    if (status) {
      boost::throw_exception(std::invalid_argument("Invalid command line option(s)"));
    } else {
      std::exit(0);
    }
  }

  void parse(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'l' :
	  license = true;
          break;
        case 'h' :
	  usage(0, std::cout);
          break;
	default :
	  usage(1);
	  break;
	}
	break;
	
      default :
	usage(1);
	break;
      }
    }
  }
};


int main(int argc, char *argv[])
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  options opts(argc, argv);
  looper::print_copyright(std::cout);
  alps::print_copyright(std::cout);
  if (opts.license) looper::print_license(std::cout);

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
