#ifndef SSE_LOOP_MATRIX_H
#define SSE_LOOP_MATRIX_H

#include <alps/parameters.h>
#include <alps/model.h>
#include <alps/lattice.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>

template <class T = double, class M = boost::numeric::ublas::matrix<T> >
class LocalMatrices
{
public:
  typedef T value_type;
  typedef M matrix_type;

  LocalMatrices(const alps::Parameters& p) :
    params_(p), graph_factory_(p), models_(p), built_(false) {}

  bool built() const { return built_; }
  void build();

  matrix_type site_matrix(int t) const { return site_matrices_[t]; }
  const matrix_type& site_matrix(int t) const { return site_matrices_[t]; }
  matrix_type bond_matrix(int t) const { return bond_matrices_[t]; }
  const matrix_type& bond_matrix(int t) const { return bond_matrices_[t]; }
  
private:
  typedef alps::graph_factory<>::graph_type graph_type;
  
  const graph_type& lattice() const { return graph_factory_.graph();}
  
  alps::Parameters parms_;
  alps::graph_factory<graph_type> graph_factory_;
  alps::ModelLibrary models_;

  bool built_;
  std::vector<matrix_type> site_matrices_;
  std::vector<matrix_type> bond_matrices_;
};

template <class T, class M>
void LocalMatrices<T,M>::build() const
{
  // get Hamilton operator
  alps::HamiltonianDescriptor<short> ham(models_.hamiltonian(parms_["MODEL"]));
  alps::Parameters p(parms_);
  p.copy_undefined(ham.default_parameters());
  ham.set_parameters(p);
  
  // get all site matrices
  alps::property_map<alps::site_type_t, graph_type, int>::const_type
  site_type(alps::get_or_default(alps::site_type_t(), lattice(), 0));

  std::map<int,boost::multi_array<T,2> > site_matrix;
  std::map<int,bool> site_visited;
  for (boost::graph_traits<graph_type>::vertex_iterator it=boost::vertices(lattice()).first; 
       it!=boost::vertices(lattice()).second ; ++it)
    if (!site_visited[site_type[*it]]) {
      int type=site_type[*it];
      std::cout << "Creating site matrix for type " << type << "\n";
      site_visited[type]=true;
      site_matrix.insert(std::make_pair(type,ham.site_term(type).template matrix<T>(
                ham.basis().site_basis(type),models_.simple_operators(),p)));
    }
  
  // get all bond matrices
  alps::property_map<alps::bond_type_t,  graph_type, int>::const_type
  bond_type(alps::get_or_default(alps::bond_type_t(),lattice(),0));

  std::map<int,boost::multi_array<T,4> > bond_matrix;
  std::map<boost::tuple<int,int,int>,bool> bond_visited;
  for (boost::graph_traits<graph_type>::edge_iterator it=boost::edges(lattice()).first; 
       it!=boost::edges(lattice()).second ; ++it) {
    int btype  = bond_type[*it];
    int stype1 = site_type[boost::source(*it,lattice())];
    int stype2 = site_type[boost::target(*it,lattice())];
    if (!bond_visited[boost::make_tuple(btype,stype1,stype2)]) {
      std::cout << "Creating bond matrix for type " << btype << "\n";
      bond_visited[boost::make_tuple(btype,stype1,stype2)]=true;
      bond_matrix.insert(std::make_pair(btype,ham.bond_term(btype).template matrix<T>(
                              ham.basis().site_basis(stype1),ham.basis().site_basis(stype2),
			      models_.simple_operators(),p)));
    }
  }
  
  // create basis set
  std::cout << "Creating basis set\n";
  alps::BasisStatesDescriptor<short> basis(ham.basis(),lattice());
  alps::BasisStates<short> states(basis);
  
  // build matrix
  
  std::cout << "Creating matrix\n";
  //matrix_.resize(boost::extents[states.size()][states.size()]);
  matrix_.resize(states.size(),states.size());

  // loop basis states
  for (int i=0;i<states.size();++i) {
    // loop over sites
    int s=0;
    for (boost::graph_traits<graph_type>::vertex_iterator it=boost::vertices(lattice()).first; 
      it!=boost::vertices(lattice()).second ; ++it,++s) {
      // get site basis index
      int is=basis[s].index(states[i][s]);
      // loop over target index
      std::vector<alps::StateDescriptor<short> > state=states[i];
      for (int js=0;js<basis[s].size();++js)
        if (site_matrix[site_type[*it]][is][js]) {
        // build target state
          state[s]=basis[s][js];
	// lookup target state
	  int j = states.index(state);
	// set matrix element
	  //matrix_[i][j]+=site_matrix[site_type[*it]][is][js];
	  matrix_(i,j)+=site_matrix[site_type[*it]][is][js];
	}
    }
    
    // loop over bonds
    for (boost::graph_traits<graph_type>::edge_iterator it=boost::edges(lattice()).first; 
      it!=boost::edges(lattice()).second ; ++it,++s) {
      // get site basis index
      int s1=boost::source(*it,lattice());
      int s2=boost::target(*it,lattice());
      int is1=basis[s1].index(states[i][s1]);
      int is2=basis[s2].index(states[i][s2]);
      // loop over target indices
      std::vector<alps::StateDescriptor<short> > state=states[i];
      for (int js1=0;js1<basis[s1].size();++js1)
        for (int js2=0;js2<basis[s2].size();++js2)
        if (bond_matrix[bond_type[*it]][is1][is2][js1][js2]) {
        // build target state
          state[s1]=basis[s1][js1];
          state[s2]=basis[s2][js2];
	// lookup target state
	  int j = states.index(state);
	// set matrix element
	  //matrix_[i][j]+=bond_matrix[bond_type[*it]][is1][is2][js1][js2];
	  matrix_(i,j)+=bond_matrix[bond_type[*it]][is1][is2][js1][js2];
	}
    }
  }  
  built_ = true;
}

#endif // SSE_LOOP_MATRIX_H
