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

template<class GRAPH>
double initialize(looper::random_choice<>& dist,
		  std::vector<looper::local_graph>& gtab,
		  const GRAPH& rlat,
		  const looper::virtual_lattice<GRAPH>& vlat,
		  const looper::model_parameter& mp)
{
  typedef typename alps::graph_traits<GRAPH>::site_iterator site_iterator;
  typedef typename alps::graph_traits<GRAPH>::bond_iterator bond_iterator;

  gtab.clear();
  std::vector<double> weight;
  double rho = 0;
  site_iterator si, si_end;
  for (boost::tie(si, si_end) = alps::sites(rlat); si != si_end; ++si) {
    looper::site_weight sw(mp.site(*si, rlat));
    site_iterator vsi, vsi_end;
    for (boost::tie(vsi, vsi_end) = virtual_sites(vlat, rlat, *si);
	 vsi != vsi_end; ++vsi)
      for (int g = 1; g <= 3; ++g)
	if (alps::is_nonzero<1>(sw.v[g])) {
	  gtab.push_back(looper::site_graph(boost::get(
            looper::site_index_t(), vlat.graph(), *vsi), g));
	  weight.push_back(sw.v[g]);
	  rho += sw.v[g];
	}
  }
  bond_iterator bi, bi_end;
  for (boost::tie(bi, bi_end) = alps::bonds(rlat); bi != bi_end; ++bi) {
    looper::bond_weight bw(mp.bond(*bi, rlat));
    bond_iterator vbi, vbi_end;
    for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat, rlat, *bi);
	 vbi != vbi_end; ++vbi)
      for (int g = 1; g <= 4; ++g)
	if (alps::is_nonzero<1>(bw.v[g])) {
	  gtab.push_back(looper::bond_graph(boost::get(
            looper::bond_index_t(), vlat.graph(), *vbi), g));
	  weight.push_back(bw.v[g]);
	  rho += bw.v[g];
	}
  }
  if (mp.has_d_term())
    for (boost::tie(si, si_end) = alps::sites(rlat); si != si_end; ++si) {
      looper::bond_weight bw(mp.site(*si, rlat));
      bond_iterator vbi, vbi_end;
      for (boost::tie(vbi, vbi_end) = virtual_bonds(vlat, rlat, *si);
	   vbi != vbi_end; ++vbi)
	for (int g = 1; g <= 4; ++g)
	  if (alps::is_nonzero<1>(bw.v[g])) {
	    gtab.push_back(looper::bond_graph(boost::get(
              looper::bond_index_t(), vlat.graph(), *vbi), g));
	    weight.push_back(bw.v[g]);
	    rho += bw.v[g];
	  }
    }
  dist.init(weight);

  return rho;
}

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
      mcs_sweep_(p["SWEEPS"]),
      vlat_(), rho_(), gtab_(), 
      r_graph(*super_type::engine_ptr, looper::random_choice<>()),
      mcs_(0)
  {
    //
    // setup model
    //

    looper::model_parameter mp(p, *this);
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
    nvsites_ = alps::num_sites(vlat_);
    nvbonds_ = alps::num_bonds(vlat_);

    //
    // setup graph table
    //

    rho_ = initialize(r_graph.distribution() , gtab_, rlat(), vlat(), mp);

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
  {
    return is_thermalized() ? (double(mcs_) / mcs_sweep_.min()) : 0.;
  }
  unsigned int mcs() const { return mcs_; }

  const graph_type& rlat() const { return super_type::graph(); }
  const looper::virtual_lattice<graph_type>& vlat() const { return vlat_; }
  unsigned int nvsites() const { return nvsites_; }
  unsigned int nvbonds() const { return nvbonds_; }

  double rho() const { return rho_; }
  const std::vector<looper::local_graph>& graph_table() const { return gtab_; }

  virtual void save(alps::ODump& od) const {
    super_type::save(od);
    od << mcs_;
  }
  virtual void load(alps::IDump& id) {
    super_type::load(id);
    id >> mcs_;
    if (super_type::where.empty()) super_type::measurements.compact();
  }

private:
  unsigned int mcs_therm_;
  looper::integer_range<unsigned int> mcs_sweep_;
  looper::virtual_lattice<graph_type> vlat_;
  unsigned int nvsites_, nvbonds_;
  double rho_;
  std::vector<looper::local_graph> gtab_;
  boost::variate_generator<alps::buffered_rng_base&, looper::random_choice<> >
    r_graph;

  // to be dumped/restored
  unsigned int mcs_;
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
  typedef looper::cluster_fragment cluster_fragment;
  typedef looper::cluster_info<qmc_type> qmc_type;

  qmc_worker(const alps::ProcessList& w, const alps::Parameters& p, int n)
    : super_type(w, p, n),
      operators(0), spins(nvsites, 0 /* all up */), operators_prev(),
      fragments(), current(nvsites), clusters(),
      r_time(*engine_ptr, boost::exponential_distribution<>(rho()))
  {}

  void dostep()
  {
    super_type::dostep();

    //
    // diagonal update and cluster construction
    //

    // initialize cluster information (setup cluster fragments)
    fragments.resize(0); fragments.resize(nvsites);
    for (int s = 0; s < nvsites; ++s) current[s] = s;
    
    // initialize operator information
    std::swap(operators, operators_prev); operators.resize(0);

    double t = r_time();

    for (std::vector<local_operator>::iterator opi = operators_prev.begin();
         t < beta || opi != operators_prev.end();) {

      // diagonal update
      if (opi == operators_prev.end() || t < opi->time) {
        // insert diagonal operator
	local_graph g = graph_table(r_graph());
	if (is_bond(g) && is_compatible(g, 
					t += r_time();
        } else {
          // insert diagonal bond operator if spins are anti-parallel
          unsigned int b = r_bond();
          if (spins[super_type::source(super_type::bond(b))] !=
              spins[super_type::target(super_type::bond(b))]) {
            operators.push_back(local_operator_t(bond_diagonal, b, t));
            t += r_time();
          } else {
            t += r_time();
            continue;
          }
        }
      } else {
        if (opi->is_diagonal()) {
          // remove diagonal operators with probability one (= nothing to do)
          ++opi;
          continue;
        } else {
          operators.push_back(*opi);
          ++opi;
        }
      }

      // building up clusters
      std::vector<local_operator_t>::reverse_iterator oi = operators.rbegin();
      if (oi->is_site()) {
        unsigned int s = oi->location;
        oi->lower_loop = current[s];
        oi->upper_loop = current[s] = add(fragments);
        if (oi->is_offdiagonal()) spins[s] ^= 1;
      } else {
        unsigned int b = oi->location;
        unsigned int s0 = super_type::source(super_type::bond(b));
        unsigned int s1 = super_type::target(super_type::bond(b));
        oi->lower_loop = unify(fragments, current[s0], current[s1]);
        oi->upper_loop = current[s0] = current[s1] = add(fragments);
        if (oi->is_offdiagonal()) {
          spins[s0] ^= 1;
          spins[s1] ^= 1;
        }
      }
    }

    // connect bottom and top cluster fragments
    for (int s = 0; s < nsites; ++s) unify(fragments, s, current[s]);

    //
    // cluster flip
    //

    // assign cluster id & determine if clusters are to be flipped
    clusters.resize(0);
    for (std::vector<fragment_t>::iterator ci = fragments.begin();
         ci != fragments.end(); ++ci)
      if (ci->is_root()) {
        ci->id = clusters.size();
        clusters.push_back(cluster_t(r_bit()));
      }

    // 'flip' operators & do improved measurements
    for (std::vector<local_operator_t>::iterator oi = operators.begin();
         oi != operators.end(); ++oi) {
      int id_l = root(fragments, oi->lower_loop).id;
      int id_u = root(fragments, oi->upper_loop).id;
      if (oi->is_site()) {
        unsigned int s = oi->location;
        clusters[id_l].mag += (1 - 2 * spins[s]) * oi->time;
        clusters[id_l].length += oi->time;
        if (oi->is_offdiagonal()) spins[s] ^= 1;
        clusters[id_u].mag -= (1 - 2 * spins[s]) * oi->time;
        clusters[id_u].length -= oi->time;
      } else {
        unsigned int b = oi->location;
        unsigned int s0 = super_type::source(super_type::bond(b));
        unsigned int s1 = super_type::target(super_type::bond(b));
        // clusters[id_l].mag += 0 * oi->time;
        clusters[id_l].length += 2 * oi->time;
        if (oi->is_offdiagonal()) {
          spins[s0] ^= 1;
          spins[s1] ^= 1;
        }
        // clusters[id_u].mag -= 0 * oi->time;
        clusters[id_u].length -= 2 * oi->time;
      }
      if (clusters[id_l].to_flip ^ clusters[id_u].to_flip) oi->flip();
    }

    // flip spins & do improved measurements
    for (unsigned int s = 0; s < nsites; ++s) {
      int id = root(fragments, s).id;
      clusters[id].mag0 += (1 - 2 * spins[s]);
      clusters[id].size += 1;
      clusters[id].mag += (1 - 2 * spins[s]) * beta;
      clusters[id].length += beta;
      if (clusters[id].to_flip) spins[s] ^= 1;
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
      for (std::vector<cluster_t>::const_iterator
             pi = clusters.begin(); pi != clusters.end(); ++pi) {
        z2 += sqr(pi->mag0);
        s2 += sqr(pi->size);
        m2 += sqr(pi->mag);
        l2 += sqr(pi->length);
      }

      measurements["Energy"]
        << (- 0.25 * model_parameter.bond().jz() * nbonds
            + 0.5  * model_parameter.site().hx() * nsites
            - (double)operators.size() / beta) / nsites;
      measurements["Uniform Magnetization^2"] << 0.25 * z2 / nsites;
      measurements["Staggered Magnetization^2"] << 0.25 * s2 / nsites;
      measurements["Uniform Susceptibility"] << 0.25 * m2 / (beta * nsites);
      measurements["Staggered Susceptibility"] << 0.25 * l2 / (beta * nsites);
    }

  }

    
  }
    
  void save(alps::ODump& od) const {
    super_type::save(od);
    od << operators << spins;
  }
  void load(alps::IDump& id) {
    super_type::load(id);
    id >> operators >> spins;
  }

private:
  std::vector<local_operator> operators;
  std::vector<int> spins;

  std::vector<local_operator> operators_prev; // no checkpointing
  std::vector<cluster_fragment> fragments; // no checkpointing
  std::vector<unsigned int> current; // no checkpointing
  std::vector<cluster_info<qmc_type> > clusters; // no checkpointing

  // RNGs
  boost::variate_generator<alps::buffered_rng_base&,
                           boost::exponential_distribution<> > r_time;
};

template<class MCRUN>
class qmc_worker<looper::sse, MCRUN> : public qmc_worker_base<MCRUN>
{
public:
  typedef looper::sse qmc_type;
  typedef qmc_worker_base<MCRUN> super_type;
  typedef looper::local_operator<qmc_type>  local_operator;

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
