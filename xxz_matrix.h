#ifndef LOOPER_XXZ_MATRIX_H
#define LOOPER_XXZ_MATRIX_H

#include <alps/parameters.h>
#include <alps/model.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace looper {

template <class T = double, class M = boost::numeric::ublas::matrix<T> >
class xxz_matrix
{
public:
  typedef T value_type;
  typedef M matrix_type;

  template<class I>
  xxz_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     double e0, double jxy, double jz) : matrix_()
  { build(s0, s1, e0, jxy, jz); }
  
  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  template<class I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     double e0, double jz, double jxy)
  {
    typedef alps::half_integer<I> half_integer_type;
    const half_integer_type u(1);
    
    // set matrix dimension
    int dim = (s0.get_twice() + 1) * (s1.get_twice() + 1);
    matrix_.resize(dim, dim);

    // diagonal elements: jz sz0 sz1
    for (half_integer_type sz0 = s0; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1; --sz1) {
	matrix_[index(s0, s1, sz0, sz1)][index(s0, s1, sz0, sz1)] = 
	  e0 - jz * double(sz0) * double(sz1);
      }
    }

    // off-diagonal elements: jxy s0+ s1- / 2
    for (half_integer_type sz0 = s0 - 1; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1 + 1; --sz1) {
	matrix_[index(s0, s1, sz0 + 1, sz1 - 1)][index(s0, s1, sz0, sz1)] =
	  - 0.5 * jxy *
	  std::sqrt(double(s0 - sz0) * double(s0 + sz0 + 1)) * 
	  std::sqrt(double(s1 + sz1) * double(s1 - sz1 + 1));
      }
    }

    // off-diagonal elements: jxy s0- s1+ / 2
    for (half_integer_type sz0 = s0; sz0 >= -s0 + 1; --sz0) {
      for (half_integer_type sz1 = s1 - 1; sz1 >= -s1; --sz1) {
	matrix_[index(s0, s1, sz0 - 1, sz1 + 1)][index(s0, s1, sz0, sz1)] =
	  - 0.5 * jxy *
	  std::sqrt(double(s0 + sz0) * double(s0 - sz0 + 1)) * 
	  std::sqrt(double(s1 - sz1) * double(s1 + sz1 + 1));
      }
    }
  }
  
protected:
  xxz_matrix();

  template<class I>
  I index(const alps::half_integer<I>& s0,
	  const alps::half_integer<I>& s1,
	  const alps::half_integer<I>& sz0,
	  const alps::half_integer<I>& sz1) const
  {
    return s0.distance(sz0) * (s1.get_twice() + 1) + s1.distance(sz1);
  }

private:
  matrix_type matrix_;
};

template <class I, class M>
inline boost::tuple<typename M::value_type, M::value_type, M::value_type>
fit_xxz_matrix(const alps::half_integer<I>& s0, 
	       const alps::half_integer<I>& s1, 
	       const M& m)
{
  return boost::make_tuple(0, 0, 0);
}

} // end namespace looper

#endif // LOOPER_XXZ_MATRIX_H
