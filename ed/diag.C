#include <alps/parameterlist.h>
#include <alps/lattice.h>
#include <alps/model.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef boost::numeric::ublas::vector<double> vector_type;
  typedef boost::numeric::ublas::matrix<double> matrix_type;
  typedef boost::numeric::ublas::vector<double> diagonal_matrix_type;

  alps::ParameterList parameterlist;
  std::cin >> parameterlist;
  for (alps::ParameterList::const_iterator p = parameterlist.begin();
       p != parameterlist.end(); ++p) {
    std::cout << *p;

    alps::graph_helper<> lattice(*p);
    alps::model_helper<> model(*p);

    alps::BasisStates<short>
      basis_set(alps::BasisStatesDescriptor<short>(model.basis(),
						   lattice.graph()));
    int dim = basis_set.size();
    std::cout << "dimension   = " << dim << std::endl;
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}
