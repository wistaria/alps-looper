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

#ifndef LOOPER_LOCATION_H
#define LOOPER_LOCATION_H

#include <alps/osiris.h>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class location
{
public:
  location(int pos = 0, bool is_bond = true)
    : loc_(pos << 1 | (is_bond ? 1 : 0)) {}
  int pos() const { return loc_ >> 1; }
  bool is_bond() const { return loc_ & 1; }
  bool is_site() const { return !is_bond(); }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
  static location bond_location(int pos) { return location(pos, true); }
  static location site_location(int pos) { return location(pos, false); }
private:
  int loc_;
};

inline int pos(const location& loc) { return loc.pos(); }
inline bool is_bond(const location& loc) { return loc.is_bond(); }
inline bool is_site(const location& loc) { return loc.is_site(); }

// optimized version for models with bond terms only

class location_bond
{
public:
  location_bond(int pos = 0, bool is_bond = true) : loc_(pos)
  {
    if (!is_bond)
      boost::throw_exception(std::invalid_argument("location_bond"));
  }
  int pos() const { return loc_; }
  static bool is_bond() { return true; }
  static bool is_site() { return false; }
  static location_bond bond_location(int pos) { return location_bond(pos); }
  static location_bond site_location(int)
  {
    boost::throw_exception(std::logic_error("location_bond"));
    return location_bond();
  }
  void save(alps::ODump& dump) const { dump << loc_; }
  void load(alps::IDump& dump) { dump >> loc_; }
private:
  int loc_;
};

inline int pos(const location_bond& loc) { return loc.pos(); }
inline bool is_bond(const location_bond&) { return true; }
inline bool is_site(const location_bond&) { return false; }

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dump, const looper::location& loc)
{ loc.save(dump); return dump; }

inline alps::IDump& operator>>(alps::IDump& dump, looper::location& loc)
{ loc.load(dump); return dump; }

inline alps::ODump& operator<<(alps::ODump& dump,
                               const looper::location_bond& loc)
{ loc.save(dump); return dump; }

inline alps::IDump& operator>>(alps::IDump& dump, looper::location_bond& loc)
{ loc.load(dump); return dump; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_LOCATION_H
