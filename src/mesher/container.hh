#ifndef CONTAINER_HH
#define CONTAINER_HH

#include <cassert>
#include <vector>
#include <algorithm>

template <class T, class Allocator = std::allocator<T>>
class unique_vector
{
protected:

  typedef std::vector<T, Allocator> vector_type;
  typedef typename vector_type::iterator iterator_type;
  typedef typename vector_type::const_iterator const_iterator_type;

public:

  inline const size_t size() const
  { return v_.size(); }

  inline const T &at(const size_t _i) const
  { assert(_i < v_.size()); return v_[_i]; }

  inline T &at(const size_t _i)
  { assert(_i < v_.size()); return v_[_i]; }

  inline const T &operator[](const size_t _i) const
  { return at(_i); }

  inline T &operator[](const size_t _i)
  { return at(_i); }

  inline iterator_type begin()
  { return v_.begin(); }

  inline const_iterator_type begin() const
  { return v_.begin(); }

  inline iterator_type end()
  { return v_.end(); }

  inline const_iterator_type end() const
  { return v_.end(); }

  inline bool exist(const T &_val) const
  { return std::find(v_.begin(), v_.end(), _val) != v_.end(); }

  inline void push_back(const T &_val)
  { if (!exist(_val)) v_.push_back(_val); }

  inline void pop_back()
  { assert(!v_.empty()); v_.pop_back(); }

  inline void clear()
  { v_.clear(); }

protected:

  // unique entries
  std::vector<T, Allocator> v_;
};

#endif // CONTAINER_HH