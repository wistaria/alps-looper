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

#include "alps/scheduler/montecarlo.h"
#include <vector>
#include <algorithm>
#include <valarray>

namespace looper {

template<typename T>
class histogram
{
public:
  histogram() : data_(), left_(0), right_(0) {}

  void resize(unsigned int newsize, unsigned int newleft = 0)
  {
    data_.resize(newsize);
    left_ = newleft;
    right_ = left_ + size() - 1;
    fill(0);
  }

  void fill(T const& x) { std::fill(data_.begin(), data_.end(), x); }

  void subtract()
  {
    T x = data_[0];
    for (unsigned int i = 0; i < data_.size(); ++i) data_[i] -= x;
  }

  std::valarray<T> getvalarray(unsigned int min, unsigned int max) const
  {
    std::valarray<T> dval;
    dval.resize(max - min + 1);
    int count = 0;
    for (int i = min-left(); i <= max-left(); ++i, ++count)
      dval[count] = data_[i];
    return dval;
  }

  const T& operator[](unsigned int i) const { return data_[i - left()]; }
  T& operator[](unsigned int i) { return data_[i - left()]; }

  unsigned int left() const { return left_; }
  unsigned int right() const { return right_; }
  unsigned int size() const { return data_.size(); }
  double flatness() const
  {
    double av = std::accumulate(data_.begin(), data_.end(), 0.);
    av /= size();
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
  operator<<(T const& x) { hist_[pos_] += x; }
private:
  histogram_t& hist_;
  unsigned int pos_;
};

template<typename T>
class histogram_set;

} // end namespace loper

#endif // LOOPER_HISTOGRAM
