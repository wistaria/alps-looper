#ifndef LOOPER_XXZ_MATRIX_H
#define LOOPER_XXZ_MATRIX_H

#include <alps/parameters.h>
#include <alps/model.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace looper {

template <class T = double, class M = boost::numeric::ublas::matrix<T> >
class xxz_matrix
{
public:
  typedef T value_type;
  typedef M matrix_type;
  typedef typename matrix_type::size_type size_type;

  xxz_matrix() : matrix_() {}
  xxz_matrix(const xxz_matrix& m) : matrix_(m.matrix_) {}
  template<class I>
  xxz_matrix(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     double e0, double jxy, double jz) : matrix_()
  { build(s0, s1, e0, jxy, jz); }
  
  // access to matrix
  matrix_type& matrix() { return matrix_; }
  const matrix_type& matrix() const { return matrix_; }

  // access to to rows
  typename matrix_type::matrix_row_type operator[](size_type i) {
    return matrix_[i];
  }
  typename matrix_type::const_matrix_row_type operator[](size_type i) const {
    return matrix_[i];
  }

  template<class I>
  void build(const alps::half_integer<I>& s0, const alps::half_integer<I>& s1,
	     double e0, double jxy, double jz)
  {
    typedef alps::half_integer<I> half_integer_type;
    
    // set matrix dimension
    int dim = (s0.get_twice()+1) * (s1.get_twice()+1);
    matrix_.resize(dim, dim);
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j)
	matrix_[i][j] = value_type(0);

    // diagonal elements: jz sz0 sz1
    for (half_integer_type sz0 = s0; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1; --sz1) {
	matrix_[index(s0, s1, sz0, sz1)][index(s0, s1, sz0, sz1)] = 
	  e0 - jz * double(sz0) * double(sz1);
      }
    }

    // off-diagonal elements: jxy s0+ s1- / 2
    for (half_integer_type sz0 = s0-1; sz0 >= -s0; --sz0) {
      for (half_integer_type sz1 = s1; sz1 >= -s1+1; --sz1) {
	matrix_[index(s0, s1, sz0+1, sz1-1)][index(s0, s1, sz0, sz1)] =
	  - 0.5 * jxy *
	  std::sqrt(double(s0-sz0) * double(s0+sz0+1)) * 
	  std::sqrt(double(s1+sz1) * double(s1-sz1+1));
      }
    }

    // off-diagonal elements: jxy s0- s1+ / 2
    for (half_integer_type sz0 = s0; sz0 >= -s0+1; --sz0) {
      for (half_integer_type sz1 = s1-1; sz1 >= -s1; --sz1) {
	matrix_[index(s0, s1, sz0-1, sz1+1)][index(s0, s1, sz0, sz1)] =
	  - 0.5 * jxy *
	  std::sqrt(double(s0+sz0) * double(s0-sz0+1)) * 
	  std::sqrt(double(s1-sz1) * double(s1+sz1+1));
      }
    }
  }
  
protected:
  template<class I>
  I index(const alps::half_integer<I>& s0,
	  const alps::half_integer<I>& s1,
	  const alps::half_integer<I>& sz0,
	  const alps::half_integer<I>& sz1) const
  {
    return s0.distance(sz0) * (s1.get_twice()+1) + s1.distance(sz1);
  }

private:
  matrix_type matrix_;
};

template <class T, class M>
inline std::ostream& operator<<(std::ostream& os, const xxz_matrix<T, M>& m)
{
  os << m.matrix();
  return os;
}

//
// fit a matrix to xxz_matrix
//

template <class I, class M>
inline boost::tuple<bool, typename M::value_type, M::value_type, M::value_type>
fit_xxz_matrix(const alps::half_integer<I>& s0, 
	       const alps::half_integer<I>& s1, 
	       const M& mat, typename M::value_type tol = 1.0e-10)
{
  typedef M matrix_type;
  typedef typename M::value_type value_type;
  typedef boost::numeric::ublas::scalar_vector<value_type> vector_type;

  int dim = (s0.get_twice()+1) * (s1.get_twice()+1);
  xxz_matrix<value_type> m1(s0, s1, 0, 1, 1);

  // e0
  value_type e0 = 0;
  for (int i = 0; i < dim; ++i) e0 += mat[i][i];
  e0 /= dim;

  // jz
  value_type jz = (mat[0][0] - e0) / m1[0][0];

  // jxy
  double jxy = 0;
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if ((i != j) && (m1[i][j] != 0)) {
	jxy = mat[i][j] / m1[i][j];
	break;
      }
    }
    if (jxy != 0) break;
  }

  // check
  bool success = true;
  xxz_matrix<value_type> m(s0, s1, e0, jxy, jz);
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      if (std::abs(mat[i][j] - m[i][j]) > tol) success = false;
    }
  }

  return boost::make_tuple(success, e0, jxy, jz);
}

} // end namespace looper

#endif // LOOPER_XXZ_MATRIX_H
