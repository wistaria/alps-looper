/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
* 
* Copyright (C) 1997-2003 by Synge Todo <wistaria@comp-phys.org>
* 
* This software is published under the ALPS Application License; you can use,
* redistribute and/or modify this software under the terms of the license,
* either version 1 or (at your option) any later version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Library; see the file LICENSE. If not, the license is also
* available from http://alps.comp-phys.org/.
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

// $Id: exact_diagonalization.h 557 2003-11-12 12:33:19Z wistaria $

#ifndef LOOPER_EXACT_DIAGONALIZATION_H
#define LOOPER_EXACT_DIAGONALIZATION_H

#include <looper/lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace looper {

template<class G, class M, class T = double>
struct exact_diagonalization
{
  typedef G graph_type;
  typedef M model_type;
  typedef T value_type;
  typedef boost::numeric::ublas::vector<value_type> vector_type;
  typedef boost::numeric::ublas::matrix<value_type> matrix_type;

  struct parameter_type
  {


  };

  struct config_type
  {

  };
};

} // end namespace looper

#endif // LOOPER_EXACT_DIAGONALIZATION_H
