#pragma once

#include <cassert>
#include <utility>

#include "global.hpp"

namespace eh
{
namespace rtree
{

/*
 * static_vector is a fixed-size vector container.
 * T element type
 * N number of elements
 * No Copy or Move Assignment defined since this would only be used in
 * child-node container
 */
template <typename T, size_type N>
class static_vector
{
public:
  using value_type = T;
  using size_type = ::eh::rtree::size_type;

protected:
  alignas(T) char data_[sizeof(T) * N];
  size_type size_;

public:
  static_vector()
      : size_(0)
  {
  }
  static_vector(static_vector const& rhs)
      : size_(0)
  {
    for (size_type i = 0; i < rhs.size(); ++i)
    {
      push_back(rhs[i]);
    }
  }
  // element-wise move construct.
  // Note that rhs.size() is not changed.
  static_vector(static_vector&& rhs)
      : size_(0)
  {
    for (size_type i = 0; i < rhs.size(); ++i)
    {
      push_back(std::move(rhs[i]));
    }
  }
  static_vector& operator=(static_vector const& rhs)
  {
    clear();
    for (size_type i = 0; i < rhs.size(); ++i)
    {
      push_back(rhs[i]);
    }
    return *this;
  }
  // element-wise move assignment
  // Note that rhs.size() is not changed.
  static_vector& operator=(static_vector&& rhs)
  {
    clear();
    for (size_type i = 0; i < rhs.size(); ++i)
    {
      push_back(std::move(rhs[i]));
    }
    return *this;
  }

  ~static_vector()
  {
    clear();
  }

  value_type& at(size_type i)
  {
    return reinterpret_cast<value_type*>(data_)[i];
  }
  value_type const& at(size_type i) const
  {
    return reinterpret_cast<value_type const*>(data_)[i];
  }
  value_type& operator[](size_type i)
  {
    return at(i);
  }
  value_type const& operator[](size_type i) const
  {
    return at(i);
  }
  void clear()
  {
    for (size_type i = 0; i < size_; ++i)
    {
      at(i).~value_type();
    }
    size_ = 0;
  }

  value_type* data()
  {
    return reinterpret_cast<value_type*>(data_);
  }
  value_type const* data() const
  {
    return reinterpret_cast<value_type const*>(data_);
  }
  size_type size() const
  {
    return size_;
  }
  size_type capacity() const
  {
    return N;
  }
  bool empty() const
  {
    return size_ == 0;
  }

  value_type& front()
  {
    assert(size_ > 0);
    return at(0);
  }
  value_type const& front() const
  {
    assert(size_ > 0);
    return at(0);
  }
  value_type& back()
  {
    assert(size_ > 0);
    return at(size_ - 1);
  }
  value_type const& back() const
  {
    assert(size_ > 0);
    return at(size_ - 1);
  }
  void push_back(value_type value)
  {
    assert(size_ < N);
    new (&at(size_)) value_type(std::move(value));
    ++size_;
  }
  template <typename... Ts>
  void emplace_back(Ts&&... args)
  {
    assert(size_ < N);
    new (&at(size_)) value_type(std::forward<Ts>(args)...);
    ++size_;
  }
  void pop_back()
  {
    assert(size_ > 0);
    back().~value_type();
    --size_;
  }

  value_type* begin()
  {
    return data();
  }
  value_type const* begin() const
  {
    return data();
  }
  value_type* end()
  {
    return data() + size();
  }
  value_type const* end() const
  {
    return data() + size();
  }
};

}
}