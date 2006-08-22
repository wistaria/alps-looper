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

#ifndef LOOPER_CLUSTER_H
#define LOOPER_CLUSTER_H

#include <boost/mpl/bool.hpp>

namespace looper {

struct cluster_info
{
  cluster_info(bool t = false) : to_flip(t), weight(0), sign(0) {}
  bool to_flip;
  double weight;
  int sign;

  template<typename F, typename FIELD>
  struct field_accumulator
  {
    // typedef F cluster_fragment_t;
    template<typename D0, typename D1, typename D2>
    field_accumulator(D0 const&, D1 const&, D2 const&) {}
    void start(int, double, int, int) const {}
    void term(int, double, int, int) const {}
  };
  template<typename F>
  struct field_accumulator<F, boost::mpl::true_>
  {
    typedef F cluster_fragment_t;
    field_accumulator(std::vector<cluster_info>& cl,
                      std::vector<cluster_fragment_t> const& fr,
                      std::vector<double> const& fd)
      : clusters(cl), fragments(fr), field(fd) {}
    void start(int p, double t, int s, int c)
    { clusters[fragments[p].id].weight -= - field[s] * (0.5-c) * t; }
    void term(int p, double t, int s, int c)
    { clusters[fragments[p].id].weight += - field[s] * (0.5-c) * t; }
    std::vector<cluster_info>& clusters;
    std::vector<cluster_fragment_t> const& fragments;
    std::vector<double> const& field;
  };

  template<typename F, typename SIGN, typename IMPROVE>
  struct sign_accumulator
  {
    typedef F cluster_fragment_t;
    template<typename D0, typename D1, typename D2, typename D3>
    sign_accumulator(D0 const&, D1 const&, D2 const&, D3 const&) {}
    void bond_sign(int, int) const {}
    void site_sign(int, int) const {}
  };
  template<typename F>
  struct sign_accumulator<F, boost::mpl::true_, boost::mpl::true_>
  {
    typedef F cluster_fragment_t;
    sign_accumulator(std::vector<cluster_info>& cl,
                     std::vector<cluster_fragment_t> const& fr,
                     std::vector<int> const& bs,
                     std::vector<int> const& ss)
      : clusters(cl), fragments(fr), bond_s(bs), site_s(ss) {}
    void bond_sign(int p, int b)
    { clusters[fragments[p].id].sign += bond_s[b]; }
    void site_sign(int p, int s)
    { clusters[fragments[p].id].sign += site_s[s]; }
    std::vector<cluster_info>& clusters;
    std::vector<cluster_fragment_t> const& fragments;
    std::vector<int> const& bond_s;
    std::vector<int> const& site_s;
  };

  template<typename F, typename FIELD, typename SIGN, typename IMPROVE>
  struct accumulator
    : private field_accumulator<F, FIELD>,
      private sign_accumulator<F, SIGN, IMPROVE>
  {
    typedef F                                  cluster_fragment_t;
    typedef field_accumulator<F, FIELD>        field_accumulator_t;
    typedef sign_accumulator<F, SIGN, IMPROVE> sign_accumulator_t;

    template<typename CL, typename FR, typename FD, typename BS, typename SS>
    accumulator(CL& cl, FR const& fr, FD const& fd, BS const& bs, SS const& ss)
      : field_accumulator_t(cl, fr, fd), sign_accumulator_t(cl, fr, bs, ss) {}
    void start_s(int p, double t, int s, int c)
    {
      field_accumulator_t::start(p, t, s, c);
      sign_accumulator_t::site_sign(p, s);
    }
    void start_b(int p0, int p1, double t, int b, int s0, int s1,
                 int c0, int c1)
    {
      field_accumulator_t::start(p0, t, s0, c0);
      field_accumulator_t::start(p1, t, s1, c1);
      sign_accumulator_t::bond_sign(p0, b);
    }
    void term_s(int p, double t, int s, int c)
    {
      field_accumulator_t::term(p, t, s, c);
      sign_accumulator_t::site_sign(p, s);
    }
    void term_b(int p0, int p1, double t, int b, int s0, int s1,
                int c0, int c1)
    {
      field_accumulator_t::term(p0, t, s0, c0);
      field_accumulator_t::term(p1, t, s1, c1);
      sign_accumulator_t::bond_sign(p0, b);
    }
    void at_bot(int p, double t, int s, int c)
    { field_accumulator_t::start(p, t, s, c); }
    void at_top(int p, double t, int s, int c)
    { field_accumulator_t::term(p, t, s, c); }
  };
};

} // end namespace looper

#endif // LOOPER_CLUSTER_H
