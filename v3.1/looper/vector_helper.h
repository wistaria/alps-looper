/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2004 by Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

#ifndef LOOPER_VECTOR_HELPER_H
#define LOOPER_VECTOR_HELPER_H

#include <algorithm>
#include <deque>
#include <functional>
#include <utility>
#include <vector>

namespace looper {

namespace { struct nulltype {}; }

template<class C = nulltype>
class index_helper
{
public:
  template<class D>
  static typename D::size_type index(const D& a,
                                     const typename D::value_type * ptr)
  {
    index_helper<D> helper(a);
    return helper.index(ptr);
  }
};


//
// for std::vector
//

template<class T, class Alloc>
class index_helper<std::vector<T, Alloc> >
{
public:
  typedef std::vector<T, Alloc>           array_type;
  typedef typename array_type::value_type value_type;
  typedef typename array_type::size_type  size_type;

  index_helper(const array_type& a) : begin_() { init(a); }

  void init(const array_type& a) { if (a.size()) begin_ = &a[0]; }

  size_type index(const value_type* ptr) const { return (ptr - begin_); }

private:
  const value_type * begin_;
};


//
// for std::deque
//

template<class T, class Alloc>
class index_helper<std::deque<T, Alloc> >
{
public:
  typedef std::deque<T, Alloc>                                   array_type;
  typedef typename array_type::value_type                        value_type;
  typedef typename array_type::size_type                         size_type;
  typedef std::vector<std::pair<size_type, size_type> > map_type;

  index_helper(const array_type& a) { init(a); }

  void init(const array_type& a)
  {
    map_.clear();

    if (a.size()) {
      // start address of first chunk
      map_.push_back(std::make_pair((size_type)(&a[0]), 0));

      // find boundary of chunks
      for (int i = 1; i < a.size(); ++i)
        if ((size_type)(&a[i]) - (size_type)(&a[i-1]) != sizeof(T))
          map_.push_back(std::make_pair((size_type)(&a[i]), i));

      std::sort(map_.begin(), map_.end());
    }
  }

  size_type index(const void* ptr) const
  {
    typename map_type::const_iterator p =
      std::lower_bound(map_.begin(), map_.end(),
                       std::make_pair((size_type)ptr, 0),
                       std::less<std::pair<size_type, size_type> >());
    if (p == map_.end() || (size_type)ptr != p->first) --p;
    return p->second + ((size_type)ptr - (p->first)) / sizeof(T);
  }

private:
  map_type map_;
};


//
// alternative implementation based on pointers for deque
// (incorrect outputs for ICC??)
//

// template<class T, class Alloc>
// class index_helper<std::deque<T, Alloc> >
// {
// public:
//   typedef std::deque<T, Alloc>                                   array_type;
//   typedef typename array_type::value_type                        value_type;
//   typedef typename array_type::size_type                         size_type;
//   typedef std::vector<std::pair<const value_type *, size_type> > map_type;
//
//   index_helper(const array_type& a) { init(a); }
//
//   void init(const array_type& a)
//   {
//     map_.clear();
//
//     if (a.size()) {
//       // start address of first chunk
//       map_.push_back(std::make_pair(&a[0], 0));
//
//       // find boundary of chunks
//       for (int i = 1; i < a.size(); ++i)
//         if ((&a[i]) - (&a[i-1]) != 1) map_.push_back(std::make_pair(&a[i], i));
//
//       std::sort(map_.begin(), map_.end());
//     }
//   }
//
//   size_type index(const value_type* ptr) const
//   {
//     typename map_type::const_iterator p =
//       std::lower_bound(map_.begin(), map_.end(), std::make_pair(ptr, 0),
//                        std::less<std::pair<const value_type *, size_type> >());
//     if (p == map_.end() || ptr != p->first) --p;
//     return p->second + (ptr - (p->first));
//   }
//
// private:
//   map_type map_;
// };


//
// naive implementation for deque (but SLOW)
//

// template<class T, class Alloc>
// class index_helper<std::deque<T, Alloc> >
// {
// public:
//   typedef std::deque<T, Alloc>                                   array_type;
//   typedef typename array_type::value_type                        value_type;
//   typedef typename array_type::size_type                         size_type;
//   typedef std::vector<std::pair<const value_type *, size_type> > map_type;
//
//   index_helper(const array_type& a) { init(a); }
//
//   void init(const array_type& a) { array_ptr_ = &a; }
//
//   size_type index(const value_type* ptr) const
//   {
//     int i = 0;
//     for (typename array_type::const_iterator itr = array_ptr_->begin();
//          itr != array_ptr_->end(); ++itr) {
//       if (&*itr == ptr) break;
//       ++i;
//     }
//     return i;
//   }
//
// private:
//   const array_type * array_ptr_;
// };

} // end namespace looper

#endif // LOOPER_VECTOR_HELPER_H
