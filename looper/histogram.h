/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 2004-2008 by Stefan Wessel <wessel@comp-phys.org>,
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
#include <alps/parapack/integer_range.h>
#include <cmath>
#include <map>
#include <string>
#include <valarray>

namespace looper {

class wl_histogram {
public:
  wl_histogram() : offset_(0) {}
  template<class T>
  explicit wl_histogram(alps::integer_range<T> const& r) : offset_(r.min())
  {
    logg_.resize(r.size(), 0);
    hist_.resize(r.size(), 0);
  }

  void subtract() { double g = logg_[0]; logg_ -= g; }
  void clear() { hist_ = 0; }

  void visit(int pos, double logf, bool update_weight) {
    int i = pos - offset_;
    if (i >= 0 && i < logg_.size()) {
      hist_[i] += 1;
      if (update_weight) logg_[i] += logf;
    }
  }

  bool check_flatness(double thresh) const {
    if (thresh < 0) return true;
    double av = hist_.sum() / hist_.size();
    double var = thresh * av;
    for (int i = 0; i < hist_.size(); ++i)
      if (std::abs(hist_[i] - av) > var) return false;
    return true;
  }
  bool check_visit(int thresh) const { return (thresh == 0) || (hist_.min() >= thresh); }

  double accept_rate(int prev, int next) const {
    int ip = prev - offset_;
    int in = next - offset_;
    if (ip >= 0 || ip < logg_.size()) {
      return (in >= 0 && in < logg_.size()) ? exp(logg_[ip] - logg_[in]) : 0;
    } else if (ip < 0) {
      return (in >= ip  && in < logg_.size()) ? 1 : 0;
    } else {
      return (in <= ip  && in >= 0) ? 1 : 0;
    }
  }

  void save(alps::ODump& dp) const { dp << offset_ << logg_ << hist_; }
  void load(alps::IDump& dp) { dp >> offset_ >> logg_ >> hist_; }

  void store(alps::ObservableSet& m, std::string const& gname,
             std::string const& hname, bool /* multicanonical */) const {
    std::valarray<double> h = hist_;
    h *= (h.size() / h.sum());
    m[gname] << logg_;
    m[hname] << h;
  }

  void output(std::ostream& os = std::cout) const {
    for (int i = 0; i < logg_.size(); ++i)
      os << (offset_+i) << ' ' << logg_[i] << ' ' << hist_[i] << std::endl;
  }

private:
  int offset_;
  std::valarray<double> logg_;
  std::valarray<double> hist_;
};


template<typename T>
class histogram_descriptor {
public:
  typedef T                                   histogram_type;
  typedef typename histogram_type::value_type value_type;
  histogram_descriptor(value_type *p) : ptr_(p) {}
  histogram_descriptor& operator<<(value_type const& x) { *ptr_ += x; return *this; }
  template<typename U>
  histogram_descriptor& operator<<(std::valarray<U> const&) {
    // ignore vector observable
    return *this;
  }
private:
  value_type *ptr_;
};

template<typename T>
class histogram_set {
public:
  typedef std::valarray<T>                      histogram_type;
  typedef std::map<std::string, histogram_type> map_type;
  typedef histogram_descriptor<histogram_type>  descriptor_type;

  histogram_set() {}
  template<typename U>
  histogram_set(alps::integer_range<U> const& r) :
    offset_(r.min()), size_(r.size()), map_() {}
  histogram_set(histogram_set const& h) :
    offset_(h.offset_), size_(h.size_), map_(h.map_), pos_(h.pos_) {}

  template<typename U>
  void initialize(alps::integer_range<U> const& r) {
    offset_ = r.min();
    size_ = r.size();
    map_.clear();
  }

  void add_histogram(std::string const& name) { map_[name].resize(size_, 0); }
  bool has(std::string const& name) const { return map_.find(name) != map_.end(); }

  void set_position(unsigned int p) { pos_ = p; }
  descriptor_type operator[](std::string const& name) {
    return descriptor_type(&map_[name][pos_ - offset_]);
  }

  void save(alps::ODump& dp) const { dp << offset_ << size_ << map_; }
  void load(alps::IDump& dp) { dp >> offset_ >> size_ >> map_; }

private:
  int offset_;
  unsigned int size_;
  map_type map_;
  unsigned int pos_;
};

template<typename T>
void add_scalar_obs(histogram_set<T>& h, std::string const& name, bool = false) {
  if (!h.has(name)) h.add_histogram(name);
}

template<typename T>
void add_vector_obs(histogram_set<T>&, std::string const&, bool = false) {
  std::cerr << "WARNING: vector observable is not supported in quantum Wang Langau\n";
}

template<typename T>
void add_vector_obs(histogram_set<T>&, std::string const&,
  alps::RealVectorObservable::label_type const&, bool = false) {
  std::cerr << "WARNING: vector observable is not supported in quantum Wang Langau\n";
}

} // end namespace looper

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace looper {
#endif

inline alps::ODump& operator<<(alps::ODump& dp, looper::wl_histogram const& h)
{ h.save(dp); return dp; }

inline alps::IDump& operator>>(alps::IDump& dp, looper::wl_histogram& h)
{ h.load(dp); return dp; }

template<typename T>
alps::ODump& operator<<(alps::ODump& dp, looper::histogram_set<T> const& h)
{ h.save(dp); return dp; }

template<typename T>
alps::IDump& operator>>(alps::IDump& dp, looper::histogram_set<T>& h)
{ h.load(dp); return dp; }

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace looper
#endif

#endif // LOOPER_HISTOGRAM
