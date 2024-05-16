#pragma once

#include <algorithm>

#include "geometry_traits.hpp"
#include "global.hpp"
#include "point.hpp"

namespace eh
{
namespace rtree
{

// bounding box representation
template <typename PointType>
struct aabb_t
{
  PointType min_, max_;

  aabb_t(PointType const& p)
      : min_(p)
      , max_(p)
  {
  }
  aabb_t(PointType const& min__, PointType const& max__)
      : min_(min__)
      , max_(max__)
  {
  }
};

template <typename ArithmeticType>
struct geometry_traits<aabb_t<ArithmeticType>>
{
  using area_type = ArithmeticType;

  using AABB = aabb_t<ArithmeticType>;

  template <typename PointType>
  static bool is_inside(AABB const& aabb, PointType const& p)
  {
    return aabb.min_ <= p && p <= aabb.max_;
  }
  template <typename PointType>
  static bool is_inside(AABB const& aabb, aabb_t<PointType> const& aabb2)
  {
    return aabb.min_ <= aabb2.min_ && aabb2.max_ <= aabb.max_;
  }

  template <typename PointType>
  static bool is_overlap(AABB const& aabb, PointType const& p)
  {
    return is_inside(aabb, p);
  }
  template <typename PointType>
  static bool is_overlap(AABB const& aabb, aabb_t<PointType> const& aabb2)
  {
    if (aabb2.min_ > aabb.max_)
    {
      return false;
    }
    if (aabb.min_ > aabb2.max_)
    {
      return false;
    }
    return true;
  }

  static AABB merge(AABB const& aabb, ArithmeticType p)
  {
    return { std::min(aabb.min_, p), std::max(aabb.max_, p) };
  }
  static AABB merge(AABB const& aabb, AABB const& aabb2)
  {
    return { std::min(aabb.min_, aabb2.min_), std::max(aabb.max_, aabb2.max_) };
  }

  static area_type area(AABB const& aabb)
  {
    return aabb.max_ - aabb.min_;
  }

  // optional; used in quadratic split resolving conflict
  static AABB intersection(AABB const& aabb, AABB const& aabb2)
  {
    const auto ret_min = std::max(aabb.min_, aabb2.min_);
    return { ret_min, std::max(ret_min, std::min(aabb.max_, aabb2.max_)) };
  }
};

// traits for multi-dimension point
template <typename T, unsigned int Dim>
struct geometry_traits<aabb_t<point_t<T, Dim>>>
{
  using area_type = T;

  using Point = point_t<T, Dim>;
  using AABB = aabb_t<Point>;

  template <typename PointType>
  static bool is_inside(AABB const& aabb, PointType const& p)
  {
    return less_equal(aabb.min_, p) && less_equal(p, aabb.max_);
  }
  template <typename PointType>
  static bool is_inside(AABB const& aabb, aabb_t<PointType> const& aabb2)
  {
    return less_equal(aabb.min_, aabb2.min_)
           && less_equal(aabb.max_, aabb.max_);
    return aabb.min_ <= aabb2.min_ && aabb2.max_ <= aabb.max_;
  }

  template <typename PointType>
  static bool is_overlap(AABB const& aabb, PointType const& p)
  {
    return is_inside(aabb, p);
  }
  template <typename PointType>
  static bool is_overlap(AABB const& aabb, aabb_t<PointType> const& aabb2)
  {
    if (!less_than(aabb2.min_, aabb.max_))
    {
      return false;
    }
    if (!less_than(aabb.min_, aabb2.max_))
    {
      return false;
    }
    return true;
  }

  static AABB merge(AABB const& aabb, Point const& p)
  {
    return { min(aabb.min_, p), max(aabb.max_, p) };
  }
  static AABB merge(AABB const& aabb, AABB const& aabb2)
  {
    return { min(aabb.min_, aabb2.min_), max(aabb.max_, aabb2.max_) };
  }

  static area_type area(AABB const& aabb)
  {
    area_type ret = 1;
    for (unsigned int i = 0; i < Dim; ++i)
    {
      ret *= (aabb.max_[i] - aabb.min_[i]);
    }
    return ret;
  }

  // optional; used in quadratic split resolving conflict
  static AABB intersection(AABB const& aabb, AABB const& aabb2)
  {
    const auto ret_min = max(aabb.min_, aabb2.min_);
    return { ret_min, max(ret_min, min(aabb.max_, aabb2.max_)) };
  }

  static bool less_than(Point const& lhs, Point const& rhs)
  {
    for (unsigned int i = 0; i < Dim; ++i)
    {
      if (lhs[i] >= rhs[i])
      {
        return false;
      }
    }
    return true;
  }
  static bool less_equal(Point const& lhs, Point const& rhs)
  {
    for (unsigned int i = 0; i < Dim; ++i)
    {
      if (lhs[i] > rhs[i])
      {
        return false;
      }
    }
    return true;
  }
  static Point min(Point const& lhs, Point const& rhs)
  {
    Point ret;
    for (unsigned int i = 0; i < Dim; ++i)
    {
      ret[i] = std::min(lhs[i], rhs[i]);
    }
    return ret;
  }

  static Point max(Point const& lhs, Point const& rhs)
  {
    Point ret;
    for (unsigned int i = 0; i < Dim; ++i)
    {
      ret[i] = std::max(lhs[i], rhs[i]);
    }
    return ret;
  }
};

}
} // namespace eh, rtree