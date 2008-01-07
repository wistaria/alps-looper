/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

// default & command line options

#include <cstdlib> // for std::exit
#include <iostream>
#include <boost/lexical_cast.hpp>

struct options
{
  unsigned int length;
  double       gamma;
  double       temperature;
  unsigned int sweeps;
  unsigned int therm;
  bool has_length;
  bool has_gamma;

  options(int argc, char *argv[], bool hl, bool hg)
    // default parameters
    : length(8), gamma(3.0), temperature(0.2),
      sweeps(1 << 16), therm(sweeps >> 3),
      has_length(hl), has_gamma(hg)
  {
    for (int i = 1; i < argc; ++i) {
      switch (argv[i][0]) {
      case '-' :
        switch (argv[i][1]) {
        case 'l' :
          if (++i == argc || !has_length) usage(1);
          length = boost::lexical_cast<unsigned int>(argv[i]); break;
        case 'g' :
          if (++i == argc || !has_gamma) usage(1);
          gamma = boost::lexical_cast<double>(argv[i]); break;
        case 't' :
          if (++i == argc) usage(1);
          temperature = boost::lexical_cast<double>(argv[i]); break;
        case 'n' :
          if (++i == argc) usage(1);
          sweeps = boost::lexical_cast<unsigned int>(argv[i]);
          therm = sweeps >> 3; break;
        case 'h' :
          usage(0, std::cout); break;
        default :
          usage(1); break;
        }
        break;
      default :
        usage(1); break;
      }
    }

    if (length % 2 == 1 || temperature <= 0. || sweeps == 0) {
      std::cerr << "invalid parameter\n"; usage(1);
    }

    if (has_length)
      std::cout << "System Length             = " << length << '\n';
    if (has_gamma)
      std::cout << "Gamma                     = " << gamma << '\n';
    std::cout << "Temperature               = " << temperature << '\n'
              << "MCS for Thermalization    = " << therm << '\n'
              << "MCS for Measurement       = " << sweeps << '\n';
  }

  void usage(int status, std::ostream& os = std::cerr) const {
    os << "[command line options]\n\n";
    if (has_length)
      os << "  -l int     System Length\n";
    if (has_gamma)
      os << "  -g double  Gamma\n";
    os << "  -t double  Temperature\n"
       << "  -n int     MCS for Measurement\n"
       << "  -h         this help\n\n";
    std::exit(status);
  }
};
