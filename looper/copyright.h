/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
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

// $Id: copyright.h 554 2003-11-12 02:36:24Z wistaria $

#ifndef LOOPER_COPYRIGHT_H
#define LOOPER_COPYRIGHT_H

#include <alps/copyright.h>
#include <iostream>

namespace looper {

std::ostream& print_copyright(std::ostream& os = std::cout)
{
  os << "ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems\n"
     << "  available from http://wistaria.comp-phys.org/looper/\n"
     << "  copyright (c) 1997-2003 by Synge Todo <wistaria@comp-phys.org>\n"
     << "\n";
  return os;
}

std::ostream& print_license(std::ostream& os = std::cout)
{
  os << "Please look at the file LICENSE for the license conditions.\n";
  return os;
}

} // end namespace looper

#endif // LOOPER_COPYRIGHT_H
