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

#ifndef LOOP_WORKER_H
#define LOOP_WORKER_H

#include <alps/alea.h>
#include <alps/scheduler.h>
#include <looper.h>

template<class MCRUN = alps::scheduler::LatticeModelMCRun<looper::graph_type> >
class qmc_worker_base : public MCRUN
{
public:
  typedef MCRUN                           super_type;
  typedef typename super_type::graph_type graph_type;
  typedef typename looper::graph_traits<graph_type>::site_iterator
    site_iterator;
  typedef typename looper::graph_traits<graph_type>::bond_iterator
    bond_iterator;

  qmc_worker_base(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n),
      mcs_therm_(static_cast<unsigned int>(p["THERMALIZATION"])),
      mcs_sweep_(p["SWEEPS"]), has_hz_(false),
      vlat_(), gtab_(), otab_(),
      r_graph_(*super_type::engine_ptr, looper::random_choice<>()),
      r_time_(*super_type::engine_ptr, boost::exponential_distribution<>()),
      mcs_(0)
  {
    //
    // setup model
    //

    looper::model_parameter mp(p, *this);
    has_hz_ = mp.has_longitudinal_field();
    bool is_signed = mp.is_signed();
    bool is_classically_frustrated = mp.is_classically_frustrated();
    if (is_signed)
      std::cerr << "WARNING: model has negative signs\n";
    if (is_classically_frustrated)
      std::cerr << "WARNING: model is classically frustrated\n";

    //
    // setup virtual lattice
    //

    bool is_bipartite = alps::set_parity(super_type::graph());
    vlat_.generate(rlat(), mp, mp.has_d_term());

    //
    // setup graph table and random number generators
    //

    double rho = initialize(mp);
    r_time_.distribution() = boost::exponential_distribution<>(rho);

    //
    // init measurements
    //

    using alps::RealObservable;
    using alps::make_observable;

    if (is_signed) {
      super_type::measurements
        << RealObservable("Sign");
    }

    // unimproved super_type::measurements
    super_type::measurements
      << make_observable(
           RealObservable("Energy"), is_signed)
      << make_observable(
           RealObservable("Energy Density"), is_signed)
      << make_observable(
           RealObservable("Diagonal Energy Density"), is_signed)
      << make_observable(
           RealObservable("Energy Density^2"), is_signed)
      << make_observable(
           RealObservable("beta * Energy / sqrt(N)"), is_signed)
      << make_observable(
           RealObservable("beta * Energy^2"), is_signed)
      << make_observable(
           RealObservable("Susceptibility"), is_signed);
    if (is_bipartite)
      super_type::measurements
        << make_observable(
             RealObservable("Staggered Susceptibility"), is_signed);

    // improved measurements
    super_type::measurements
      << make_observable(
           RealObservable("Magnetization^2"), is_signed)
      << make_observable(
           RealObservable("Diagonal Energy Density (improved)"),
           is_signed);
    if (is_bipartite)
      super_type::measurements
        << make_observable(
             RealObservable("Staggered Magnetization^2"),
             is_signed);
    if (!is_classically_frustrated) {
      super_type::measurements
        << RealObservable("Uniform Generalized Magnetization^2")
        << RealObservable("Uniform Generalized Susceptibility");
      if (is_bipartite)
        super_type::measurements
          << RealObservable("Staggered Generalized Magnetization^2")
          << RealObservable("Staggered Generalized Susceptibility");
    }

  }

  virtual ~qmc_worker_base() {}

  virtual void dostep() { ++mcs_; }
  bool is_thermalized() const { return mcs_ >= mcs_therm_; }
  double work_done() const
  { return is_thermalized() ? (double(mcs_) / mcs_sweep_.min()) : 0.; }
  unsigned int mcs() const { return mcs_; }

  bool has_longitudinal_field() const { return has_hz_; }

  const graph_type& rlat() const { return super_type::graph(); }
  const looper::virtual_lattice<graph_type>& vlat() const { return vlat_; }
  unsigned int vsource(unsigned int b) const
  { return boost::source(*(looper::bonds(vlat_).first + b), vlat_.graph()); }
  unsigned int vtarget(unsigned int b) const
  { return boost::target(*(looper::bonds(vlat_).first + b), vlat_.graph()); }

  double advance() const { return r_time_(); }
  const looper::local_graph& choose_graph() const { return gtab_[r_graph_()]; }
  looper::local_graph choose_graph(const looper::location& loc) const
  {
    int g = (is_site(loc) || random() < otab_[pos(loc)]) ? 0 : 2;
    return looper::local_graph(g, loc);
  }

  virtual void save(alps::ODump& od) const {
    super_type::save(od);
    od << mcs_;
  }
  virtual void load(alps::IDump& id) {
    super_type::load(id);
    id >> mcs_;
    if (super_type::where.empty()) super_type::measurements.compact();
  }

protected:
  double initialize(const looper::model_parameter& mp)
  {
    gtab_.clear();
    otab_.clear();
    std::vector<double> weight;
    double rho = 0;
    site_iterator si, si_end;
    for (boost::tie(si, si_end) = alps::sites(rlat()); si != si_end; ++si) {
      looper::site_weight sw(mp.site(*si, rlat()));
      site_iterator vsi, vsi_end;
      for (boost::tie(vsi, vsi_end) = virtual_sites(vlat_, rlat(), *si);
           vsi != vsi_end; ++vsi)
        for (int g = 0; g <= 2; ++g)
          if (alps::is_nonzero<1>(sw.v[g])) {
            gtab_.push_back(looper::site_graph(g, boost::get(
              looper::site_index_t(), vlat_.graph(), *vsi)));
            weight.push_back(sw.v[g]);
            rho += sw.v[g];
          }
    }
    bond_iterator bi, bi_end;
    for (boost::tie(bi, bi_end) = alps::bonds(rlat()); bi != bi_end; ++bi) {
      looper::bond_weight bw(mp.bond(*bi, rlat()));
      bond_iterator vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat_, rlat(), *bi);
           vbi != vbi_end; ++vbi) {
        for (int g = 0; g <= 3; ++g)
          if (alps::is_nonzero<1>(bw.v[g])) {
            gtab_.push_back(looper::bond_graph(g, boost::get(
              looper::bond_index_t(), vlat_.graph(), *vbi)));
            weight.push_back(bw.v[g]);
            rho += bw.v[g];
          }
        if (alps::is_nonzero<1>(bw.v[0] + bw.v[2]))
          otab_.push_back(bw.v[0] / (bw.v[0] + bw.v[2]));
        else
          otab_.push_back(1);
      }
    }
    if (mp.has_d_term())
      for (boost::tie(si, si_end) = alps::sites(rlat()); si != si_end; ++si) {
        looper::bond_weight bw(mp.site(*si, rlat()));
        bond_iterator vbi, vbi_end;
        for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat_, rlat(), *si);
             vbi != vbi_end; ++vbi)
          for (int g = 0; g <= 3; ++g)
            if (alps::is_nonzero<1>(bw.v[g])) {
              gtab_.push_back(looper::bond_graph(g, boost::get(
                looper::bond_index_t(), vlat_.graph(), *vbi)));
              weight.push_back(bw.v[g]);
              rho += bw.v[g];
            }
      }
    r_graph_.distribution().init(weight);
    return rho;
  }

private:
  unsigned int mcs_therm_;
  looper::integer_range<unsigned int> mcs_sweep_;
  bool has_hz_;
  looper::virtual_lattice<graph_type> vlat_;
  std::vector<looper::local_graph> gtab_; // graph probability for diagonal conf
  std::vector<double> otab_; // graph probability for offdiagonal conf
  mutable boost::variate_generator<alps::buffered_rng_base&,
                                   looper::random_choice<> > r_graph_;
  mutable boost::variate_generator<alps::buffered_rng_base&,
                                   boost::exponential_distribution<> > r_time_;

  // to be dumped/restored
  unsigned int mcs_;
};

template<class QMC> struct cluster_info;

template<>
struct cluster_info<looper::path_integral>
{
  cluster_info(bool t = false)
    : to_flip(t), mag0(0), size(0), mag(0), length(0) {}
  bool to_flip;
  int mag0;
  int size;
  double mag;
  double length;
};

template<class QMC,
  class MCRUN = alps::scheduler::LatticeModelMCRun<looper::graph_type> >
class qmc_worker;

template<class MCRUN>
class qmc_worker<looper::path_integral, MCRUN> : public qmc_worker_base<MCRUN>
{
public:
  typedef looper::path_integral qmc_type;
  typedef qmc_worker_base<MCRUN> super_type;
  typedef looper::local_operator<qmc_type>  local_operator;
  typedef looper::union_find::node_idx cluster_fragment;
  typedef cluster_info<qmc_type> cluster_info;
  typedef typename super_type::site_iterator site_iterator;

  qmc_worker(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n),
      spins(looper::num_sites(super_type::vlat()), 0 /* all up */),
      operators(0), beta(1.0 / static_cast<double>(p["T"])),
      spins_curr(looper::num_sites(super_type::vlat())), operators_prev(),
      fragments(), current(looper::num_sites(super_type::vlat())), clusters()
  {}

  void dostep()
  {
    super_type::dostep();

    const int ns = looper::num_sites(super_type::vlat());

    //
    // diagonal update and cluster construction
    //

    std::copy(spins.begin(), spins.end(), spins_curr.begin());

    // initialize cluster information (setup cluster fragments)
    fragments.resize(0); fragments.resize(ns);
    for (int s = 0; s < ns; ++s) current[s] = s;
    int ghost = super_type::has_longitudinal_field() ? add(fragments) : 0;

    // initialize operator information
    std::swap(operators, operators_prev); operators.resize(0);

    double t = super_type::advance();

    for (std::vector<local_operator>::iterator opi = operators_prev.begin();
         t < beta || opi != operators_prev.end();) {

      // diagonal update & labeling
      if (opi == operators_prev.end() || t < opi->time()) {
        // insert diagonal operator and graph if compatible
        looper::local_graph g = super_type::choose_graph();
        if ((is_bond(g) &&
             is_compatible(g, spins_curr[super_type::vsource(pos(g))],
                           spins_curr[super_type::vtarget(pos(g))])) ||
            (is_site(g) && is_compatible(g, spins_curr[pos(g)]))) {
          operators.push_back(local_operator(g, t));
          t += super_type::advance();
        } else {
          t += super_type::advance();
          continue;
        }
      } else {
        if (opi->is_diagonal()) {
          // remove diagonal operator with probability one (= nothing to do)
          ++opi;
          continue;
        } else {
          // assign graph to offdiagonal operator
          opi->assign_graph(super_type::choose_graph(opi->loc()));
          operators.push_back(*opi);
          ++opi;
        }
      }

      // building up clusters
      std::copy(spins.begin(), spins.end(), spins_curr.begin());
      std::vector<local_operator>::reverse_iterator oi = operators.rbegin();
      if (oi->is_bond()) {
        unsigned int b = oi->pos();
        unsigned int s0 = super_type::vsource(b);
        unsigned int s1 = super_type::vtarget(b);
        boost::tie(current[s0], current[s1], oi->loop0, oi->loop1) =
          reconnect(fragments, oi->graph(), current[s0], current[s1]);
        if (oi->is_offdiagonal()) {
          spins_curr[s0] ^= 1;
          spins_curr[s1] ^= 1;
        }
      } else {
        unsigned int s = oi->pos();
        boost::tie(current[s], oi->loop0, oi->loop1) =
          reconnect(fragments, oi->graph(), current[s]);
        if (oi->is_locked()) unify(fragments, ghost, current[s]);
        if (oi->is_offdiagonal()) spins_curr[s] ^= 1;
      }
    }

    // connect bottom and top cluster fragments after random permutation
    {
      std::vector<int> r, c0, c1;
      site_iterator rsi, rsi_end;
      for (boost::tie(rsi, rsi_end) = alps::sites(super_type::rlat());
           rsi != rsi_end; ++rsi) {
        site_iterator vsi, vsi_end;
        boost::tie(vsi, vsi_end) =
          virtual_sites(super_type::vlat(), super_type::rlat(), *rsi);
        int offset = *vsi;
        int s2 = *vsi_end - *vsi;
        r.resize(s2); c0.resize(s2); c1.resize(s2);
        for (int i = 0; i < s2; ++i) {
          r[i] = i;
          c0[i] = spins[offset+i];
          c1[i] = spins_curr[offset+i];
        }
        looper::restricted_random_shuffle(r.begin(), r.end(), c0.begin(),
          c0.end(), c1.begin(), c1.end(), *super_type::engine_ptr);
        for (int i = 0; i < s2; ++i)
          unify(fragments, offset+i, current[offset+r[i]]);
      }
    }

    //
    // cluster flip
    //

    // assign cluster id & determine if clusters are to be flipped
    clusters.resize(0);
    for (std::vector<cluster_fragment>::iterator ci = fragments.begin();
         ci != fragments.end(); ++ci)
      if (ci->is_root()) {
        ci->id = clusters.size();
        clusters.push_back(cluster_info(random() < 0.5));
      }
    if (super_type::has_longitudinal_field())
      clusters[cluster_id(fragments, ghost)].to_flip = false;

    // flip operators and spins & do improved measurements
    if (!super_type::has_longitudinal_field()) {
      std::copy(spins.begin(), spins.end(), spins_curr.begin());
      for (std::vector<local_operator>::iterator oi = operators.begin();
           oi != operators.end(); ++oi) {
        int id_l = root(fragments, oi->loop0).id;
        int id_u = root(fragments, oi->loop1).id;
        if (oi->is_site()) {
          unsigned int s = oi->pos();
          clusters[id_l].mag += (1 - 2 * spins[s]) * oi->time();
          clusters[id_l].length += oi->time();
          if (oi->is_offdiagonal()) spins[s] ^= 1;
          clusters[id_u].mag -= (1 - 2 * spins[s]) * oi->time();
          clusters[id_u].length -= oi->time();
        } else {
          unsigned int b = oi->pos();
          unsigned int s0 = super_type::vsource(b);
          unsigned int s1 = super_type::vtarget(b);
          // clusters[id_l].mag += 0 * oi->time;
          clusters[id_l].length += 2 * oi->time();
          if (oi->is_offdiagonal()) {
            spins[s0] ^= 1;
            spins[s1] ^= 1;
          }
          // clusters[id_u].mag -= 0 * oi->time;
          clusters[id_u].length -= 2 * oi->time();
        }
        if (clusters[id_l].to_flip ^ clusters[id_u].to_flip) oi->flip();
      }
      for (unsigned int s = 0; s < ns; ++s) {
        int id = cluster_id(fragments, s);
        clusters[id].mag0 += (1 - 2 * spins[s]);
        clusters[id].size += 1;
        clusters[id].mag += (1 - 2 * spins[s]) * beta;
        clusters[id].length += beta;
        if (clusters[id].to_flip) spins[s] ^= 1;
      }
    } else {
      for (std::vector<local_operator>::iterator oi = operators.begin();
           oi != operators.end(); ++oi)
        if (clusters[cluster_id(fragments, oi->loop0)].to_flip ^
            clusters[cluster_id(fragments, oi->loop1)].to_flip)
          oi->flip();
      for (unsigned int s = 0; s < ns; ++s)
        if (clusters[cluster_id(fragments, s)].to_flip) spins[s] ^= 1;
    }

    //
    // measurements
    //

    {
      // accumurate loop length and magnetization
      double z2 = 0;
      double s2 = 0;
      double m2 = 0;
      double l2 = 0;
      std::vector<cluster_info>::iterator x;
      for (std::vector<cluster_info>::iterator pi = clusters.begin();
           pi != clusters.end(); ++pi) {
        z2 += looper::sqr(pi->mag0);
        s2 += looper::sqr(pi->size);
        m2 += looper::sqr(pi->mag);
        l2 += looper::sqr(pi->length);
      }

      int n = alps::num_sites(super_type::rlat());
      super_type::measurements["Energy"]
        << - (double)operators.size() / beta / n;
      super_type::measurements["Uniform Magnetization^2"] << z2 / (4 * n);
      super_type::measurements["Staggered Magnetization^2"] << s2 / (4 * n);
      super_type::measurements["Uniform Susceptibility"] << m2 / (4 * beta * n);
      super_type::measurements["Staggered Susceptibility"]
        << l2 / (4 * beta * n);
    }

  }

  void save(alps::ODump& od) const {
    super_type::save(od);
    // od << operators << spins;
  }
  void load(alps::IDump& id) {
    super_type::load(id);
    // id >> operators >> spins;
  }

private:
  std::vector<int> spins;
  std::vector<local_operator> operators;

  // no checkpointing
  double beta;
  std::vector<int> spins_curr;
  std::vector<local_operator> operators_prev;
  std::vector<cluster_fragment> fragments;
  std::vector<unsigned int> current;
  std::vector<cluster_info> clusters;
};

template<class MCRUN>
class qmc_worker<looper::sse, MCRUN> : public qmc_worker_base<MCRUN>
{
public:
  typedef looper::sse qmc_type;
  typedef qmc_worker_base<MCRUN> super_type;
  typedef looper::local_operator<qmc_type> local_operator;

  qmc_worker(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n)
  {}

  void save(alps::ODump& od) const { super_type::save(od); }
  void load(alps::IDump& id) { super_type::load(id); }
};


template<class T>
inline void accumulate(const alps::ObservableSet& m_in, T& m_out)
{
  if (m_in.has("beta * Energy / sqrt(N)") && m_in.has("beta * Energy^2")) {
    alps::RealObsevaluator obse_e = m_in["beta * Energy / sqrt(N)"];
    alps::RealObsevaluator obse_e2 = m_in["beta * Energy^2"];
    alps::RealObsevaluator eval("Specific Heat");
    eval = (obse_e2 - obse_e * obse_e);
    m_out << eval;
  }
}

inline void accumulate(alps::scheduler::MCSimulation& sim)
{ accumulate(sim.get_measurements(), sim); }

#endif // LOOP_WORKER_H
