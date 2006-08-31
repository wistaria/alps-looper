/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_INTEGER_RANGE_H
#define LOOPER_INTEGER_RANGE_H

#include <alps/expression.h>
#include <alps/osiris.h>
#include <boost/call_traits.hpp>
#include <boost/spirit/core.hpp>
#include <boost/throw_exception.hpp>
#include <complex>
#include <iosfwd>
#include <limits>
#include <stdexcept>
#include <string>

namespace looper {

//
// class tempalte integer_range
//

template<class T>
class integer_range
{
public:
  typedef T value_type;
  typedef typename boost::call_traits<value_type>::param_type param_type;

  integer_range() {}
  explicit integer_range(param_type v) : mi_(v), ma_(v) {}
  explicit integer_range(param_type vmin, param_type vmax)
    : mi_(vmin), ma_(vmax) {}
  integer_range(integer_range const& r) : mi_(r.mi_), ma_(r.ma_) {}
  integer_range(std::string const& str,
                alps::Parameters const& p = alps::Parameters(),
                param_type def_mi = std::numeric_limits<value_type>::min(),
                param_type def_ma = std::numeric_limits<value_type>::max())
    : mi_(def_mi), ma_(def_ma)
  {
    using namespace boost::spirit;
    std::string mi_str, ma_str;
    bool success;
    success = parse(str.c_str(),
        ( ch_p('[')
          >> (*(anychar_p-'['-':'-']'))[assign_a(mi_str)]
          >> ':'
          >> (*(anychar_p-'['-':'-']'))[assign_a(ma_str)]
          >> ']'
        )
      | ( ch_p('[')
          >> (*(anychar_p-'['-':'-']'))[assign_a(mi_str)][assign_a(ma_str)]
          >> ']'
        )
      | (+(anychar_p-'['-':'-']'))[assign_a(mi_str)][assign_a(ma_str)]
      , space_p).full;
    if (!success)
      boost::throw_exception(std::runtime_error("integer_range: parse error"));
    if (!mi_str.empty()) {
      int m = alps::evaluate(mi_str, p);
      if (!std::numeric_limits<value_type>::is_signed && m < 0)
        boost::throw_exception(std::runtime_error(
          "integer_range: negative minimum value not allowed"));
      mi_ = m;
    }
    if (!ma_str.empty()) {
      int m = alps::evaluate(ma_str, p);
      if (!std::numeric_limits<value_type>::is_signed && m < 0)
        boost::throw_exception(std::runtime_error(
          "integer_range: negative maximum value not allowed"));
      ma_ = m;
    }
    if (mi_ > ma_)
      boost::throw_exception(std::runtime_error("integer_range: range error"));
  }
  integer_range& operator=(param_type v)
  {
    mi_ = ma_ = v;
    return *this;
  }

  void include(param_type v)
  {
    if (v < mi_) mi_ = v;
    if (v > ma_) ma_ = v;
  }

  value_type min() const { return mi_; }
  value_type max() const { return ma_; }
  value_type size() const { return ma_ - mi_ + 1; }
  bool is_included(param_type v) const { return (v >= min()) && (v <= max()); }

  void save(alps::ODump& dp) const { dp << mi_ << ma_; }
  void load(alps::IDump& dp) { dp >> mi_ >> ma_; }

private:
  value_type mi_, ma_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

template<class T>
std::ostream& operator<<(std::ostream& os, integer_range<T> const& ir)
{ os << '[' << ir.min() << ':' << ir.max() << ']'; return os; }

template<class T>
alps::ODump& operator<<(alps::ODump& dp, looper::integer_range<T> const& ir)
{ ir.save(dp); return dp; }

template<class T>
alps::IDump& operator>>(alps::IDump& dp, looper::integer_range<T>& ir)
{ ir.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_INTEGER_RANGE_H
