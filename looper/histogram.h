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

#include <alps/osiris.h>
#include <boost/call_traits.hpp>
#include <boost/lambda/lambda.hpp>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <valarray>
#include <vector>

namespace looper {

template<typename T>
class histogram
{
public:
  histogram() : left_(0) {}

  void resize(unsigned int size, unsigned int left = 0)
  {
    data_.clear();
    data_.resize(size, 0);
    left_ = left;
  }

  void fill(typename boost::call_traits<T>::param_type x)
  { std::fill(data_.begin(), data_.end(), x); }

  void subtract()
  { std::for_each(data_.begin(), data_.end(), boost::lambda::_1 -= data_[0]); }

  std::valarray<T> get_array(unsigned int min, unsigned int max) const
  {
    std::valarray<T> dval;
    dval.resize(max - min + 1);
    int count = 0;
    for (int i = min - left_; i <= max - left_; ++i, ++count)
      dval[count] = data_[i];
    return dval;
  }

  const T& operator[](unsigned int i) const { return data_[i - left_]; }
  T& operator[](unsigned int i) { return data_[i - left_]; }

  unsigned int left() const { return left_; }
  unsigned int right() const { return left_ + size() - 1; }
  unsigned int size() const { return data_.size(); }
  double flatness() const
  {
    double av = std::accumulate(data_.begin(), data_.end(), 0.) / size();
    double diff = std::abs(data_[0] - av);
    for (int i = 1; i < size(); ++i)
      diff = std::max(diff, std::abs(data_[i] - av));
    return dip(diff, av);
  }
  T min() const
  {
    T min = data_[0];
    for (int i = 1; i < size(); ++i) if (data_[i] < min) min = data_[i];
    return min;
  }
  void save(alps::ODump& dump) const { dump << data_ << left_; }
  void load(alps::IDump& dump)
  { dump >> data_ >> left_; right_ = left_ + size() - 1; }

private:
  std::vector<T> data_;
  unsigned int left_;
  unsigned int right_;
};

template<typename T>
class histogram_descriptor
{
public:
  typedef histogram<T> histogram_t;
  histogram_descriptor(histogram_t& h, unsigned int p) : hist_(h), pos_(p) {}
  histogram_descriptor& operator<<(T const& x)
  {
    hist_[pos_] += x;
    return *this;
  }
private:
  histogram_t& hist_;
  unsigned int pos_;
};

template<typename T>
class histogram_set : public std::map<std::string, histogram<T> >
{
public:
  typedef std::map<std::string, histogram<T> > super_type;
  typedef histogram<T>                         histogram_type;
  typedef histogram_descriptor<T>              descriptor_type;

  histogram_set() : super_type() {}
  histogram_set(histogram_set const& h) : super_type(h) {}
  ~histogram_set() {}

  void add_histogram(std::string const& name)
  { super_type::operator[](name) = histogram_type(); }

  void set_index(unsigned int p) { index_ = p; }

  descriptor_type& operator[](std::string const& name)
  { return descriptor_type(super_type::operator[](name), index_); }

private:
  unsigned int index_;
};

} // end namespace loper

#endif // LOOPER_HISTOGRAM
