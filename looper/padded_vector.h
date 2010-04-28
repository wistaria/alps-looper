/*****************************************************************************
*
* ALPS/looper: multi-cluster quantum Monte Carlo algorithms for spin systems
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_PADDED_VECTOR
#define LOOPER_PADDED_VECTOR

namespace looper {

template<typename CONTAINER, unsigned int PADDING_BITS>
class padded_vector
{
private:
  typedef CONTAINER container_t;
public:
  typedef std::size_t size_type;
  typedef typename CONTAINER::value_type value_type;
  typedef typename CONTAINER::reference reference;
  typedef typename CONTAINER::const_reference const_reference;
  static const int padding_bits = PADDING_BITS;
  static const int padding_size = (1 << padding_bits);
  padded_vector() : container_() {}
  explicit padded_vector(size_type n, const_reference x = value_type()) :
    container_(n << padding_bits, x) {}
  reference operator[](size_type i) { return container_[i << padding_bits]; }
  const_reference operator[](size_type i) const { return container_[i << padding_bits]; }
  void resize(size_type n, const_reference x = value_type()) {
    container_.resize(n << padding_bits, x);
  }
private:
  container_t container_;
};

} // end namespace looper

#endif // LOOPER_PADDED_VECTOR
