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

#ifndef LOOPER_TEMPERATURE_H
#define LOOPER_TEMPERATURE_H

#include <alps/expression.h>
#include <boost/next_prior.hpp>
#include <boost/throw_exception.hpp>
#include <stdexcept>

namespace looper {

class temperature
{
public:
  temperature(alps::Parameters const& p) { init(p); }

  void init(alps::Parameters const& p)
  {
    seq_.clear();
    final_ = 0;

    if (p.defined("T")) final_ = alps::evaluate("T", p);
    if (final_ <= 0)
      boost::throw_exception(std::invalid_argument(
        "non-positive temperature"));

    int s_prev = 1;
    for (int i = 0;; ++i) {
      std::string ns = boost::lexical_cast<std::string>(i);
      if (p.defined("T_DURATION_" + ns) && p.defined("T_START_" + ns)) {
        int td = static_cast<int>(alps::evaluate("T_DURATION_" + ns, p));
        double ts = alps::evaluate("T_START_" + ns, p);
        if (td <= 0 || ts <= 0) break;
        seq_.push_back(std::make_pair(s_prev, ts));
        s_prev += td;
      } else {
        break;
      }
    }
    seq_.push_back(std::make_pair(s_prev, final_));
  }

  double initial() const { return seq_.front().second; }
  double final() const { return final_; }
  int annealing_steps() const { return seq_.back().first - 1; }
  double operator()(int step) const
  {
    if (step <= 0)
      boost::throw_exception(std::invalid_argument("invalid MCS"));
    std::vector<std::pair<int, double> >::const_iterator itr =
      std::lower_bound(seq_.begin(), seq_.end(), std::make_pair(step, 0.));
    if (itr == seq_.begin())
      return itr->second;
    else if (itr != seq_.end())
      return boost::prior(itr)->second +
        (step - boost::prior(itr)->first) *
        (itr->second - boost::prior(itr)->second) /
        (itr->first - boost::prior(itr)->first);
    else
      return final_;
  }

protected:
  temperature() {}

private:
  std::vector<std::pair<int, double> > seq_;
  double final_;
};

} // end namespace looper

#endif // LOOPER_TEMPERATURE_H
