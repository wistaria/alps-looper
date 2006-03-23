/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004-2006 by Stefan Wessel <wessel@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_HISTOGRAM_H
#define LOOPER_HISTOGRAM_H

#include "util.h"

#include <alps/osiris.h>
#include <boost/call_traits.hpp>
#include <boost/lambda/lambda.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <valarray>
#include <vector>

#if 1

namespace looper {

class wl_histogram
{
public:
  wl_histogram() : offset_(0) {}
  template<class T>
  explicit wl_histogram(integer_range<T> const& r)
    : offset_(r.min()), accept_(r.size(), 1), hist_(r.size(), 0)
  { accept_.back() = 0; }

  template<class T>
  void resize(integer_range<T> const& r)
  {
    offset_ = r.min();
    accept_.clear();
    accept_.resize(r.size(), 1);
    accept_.back() = 0;
    hist_.clear();
    hist_.resize(r.size(), 0);
  }
  void clear() { std::fill(hist_.begin(), hist_.end(), 0); }

  void visit(int nop, double factor)
  {
    int i = nop - offset_;
    if (i >= 1) accept_[i-1] /= factor;
    accept_[i] *= factor;
    hist_[i] += 1;
  }

  bool check_flatness(double thresh) const
  {
    if (thresh < 0) return true;
    int av = std::accumulate(hist_.begin(), hist_.end(), 0) / hist_.size();
    int tn = static_cast<int>(thresh * av);
    for (std::vector<int>::const_iterator itr = hist_.begin();
         itr != hist_.end(); ++itr) if (std::abs(*itr - av) >= tn) return false;
    return true;
  }
  bool check_visit(int thresh) const
  {
    if (thresh == 0) return true;
    for (std::vector<int>::const_iterator itr = hist_.begin();
         itr != hist_.end(); ++itr) if (*itr < thresh) return false;
    return true;
  }

  double accept_rate(int nop) const { return accept_[nop - offset_]; }
  int operator[](int nop) const { return hist_[nop - offset_]; }

  void save(alps::ODump& dp) const { dp << offset_ << accept_ << hist_; }
  void load(alps::IDump& dp) { dp >> offset_ >> accept_ >> hist_; }

  void output_dos(std::string const& prefix) const
  {
    double logg = 0;
    for (int i = 0; i < accept_.size(); ++i) {
      std::cout << prefix << ' ' << offset_ + i << ' ' << logg << std::endl;
      logg -= log(accept_[i]);
    }
  }

private:
  int offset_;
  std::vector<double> accept_;
  std::vector<int> hist_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::wl_histogram const& h)
{ h.save(dp); return dp; }

inline alps::IDump& operator>>(alps::IDump& dp, looper::wl_histogram& h)
{ h.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

# else

namespace looper {

class wl_histogram
{
public:
  wl_histogram() : offset_(0) {}
  template<class T>
  explicit wl_histogram(integer_range<T> const& r)
    : offset_(r.min()), logg_(r.size(), 0), hist_(r.size(), 0) {}

  template<class T>
  void resize(integer_range<T> const& r)
  {
    offset_ = r.min();
    logg_.clear();
    logg_.resize(r.size(), 0);
    hist_.clear();
    hist_.resize(r.size(), 0);
  }
  void clear() { std::fill(hist_.begin(), hist_.end(), 0); }

  void visit(int nop, double logf)
  {
    int i = nop - offset_;
    logg_[i] += logf;
    hist_[i] += 1;
  }

  bool check_flatness(double thresh) const
  {
    if (thresh < 0) return true;
    int av = std::accumulate(hist_.begin(), hist_.end(), 0) / hist_.size();
    int tn = static_cast<int>(thresh * av);
    for (std::vector<int>::const_iterator itr = hist_.begin();
         itr != hist_.end(); ++itr) if (std::abs(*itr - av) >= tn) return false;
    return true;
  }
  bool check_visit(int thresh) const
  {
    if (thresh == 0) return true;
    for (std::vector<int>::const_iterator itr = hist_.begin();
         itr != hist_.end(); ++itr) if (*itr < thresh) return false;
    return true;
  }

  double accept_rate(int nop) const
  {
    int i = nop - offset_;
    return (i < logg_.size() - 1) ? exp(logg_[i] - logg_[i+1]) : 0;
  }
  int operator[](int nop) const { return hist_[nop - offset_]; }

  void save(alps::ODump& dp) const { dp << offset_ << logg_ << hist_; }
  void load(alps::IDump& dp) { dp >> offset_ >> logg_ >> hist_; }

  void output_dos(std::string const& prefix) const
  {
    for (int i = 0; i < logg_.size(); ++i)
      std::cout << prefix << ' ' << offset_ + i << ' ' << logg_[i] << std::endl;
  }

private:
  int offset_;
  std::vector<double> logg_;
  std::vector<int> hist_;
};

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::wl_histogram const& h)
{ h.save(dp); return dp; }

inline alps::IDump& operator>>(alps::IDump& dp, looper::wl_histogram& h)
{ h.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif

namespace looper {

template<typename T>
class histogram_descriptor
{
public:
  typedef T                                   histogram_type;
  typedef typename histogram_type::value_type value_type;
  histogram_descriptor(value_type *p) : ptr_(p) {}
  histogram_descriptor& operator<<(value_type const& x)
  { *ptr_ += x; return *this; }
private:
  value_type *ptr_;
};

template<typename T>
class histogram_set
{
public:
  typedef std::vector<T>                         histogram_type;
  typedef std::map<std::string, histogram_type > map_type;
  typedef histogram_descriptor<histogram_type>   descriptor_type;

  histogram_set() {}
  template<typename U>
  histogram_set(integer_range<U> const& r) :
    offset_(r.min()), size_(r.size()), map_() {}
  histogram_set(histogram_set const& h) :
    offset_(h.offset_), size_(h.size_), map_(h.map_), pos_(h.pos_) {}

  template<typename U>
  void initialize(integer_range<U> const& r)
  {
    offset_ = r.min();
    size_ = r.size();
    map_.clear();
  }

  void add_histogram(std::string const& name)
  { map_[name] = histogram_type(size_, 0); }
  bool has(std::string const& name) const
  { return map_.find(name) != map_.end(); }

  void set_position(unsigned int p) { pos_ = p; }
  descriptor_type operator[](std::string const& name)
  { return descriptor_type(&map_[name][pos_ - offset_]); }

private:
  int offset_;
  unsigned int size_;
  map_type map_;
  unsigned int pos_;
};


template<typename T>
void add_measurement(histogram_set<T>& h, std::string const& name,
                     bool /* is_signed */ = false)
{ if (!h.has(name)) h.add_histogram(name); }

} // end namespace loper

#endif // LOOPER_HISTOGRAM
