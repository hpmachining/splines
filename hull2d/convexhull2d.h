// ConvexHull2D
// Template class to calculate the 2d convex hull from a 2d point cloud.
// The ConvexHull2d class takes a std::vector<T>, where <T> is any
// class or structure that has public "x()" and "y()" methods representing coordinates.
//
// Uses Andrew's Monotone Chain algorithm
// Algorithm based on code found at this web site:
// https://en.m.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
// 
// Example:
//    std::vector<MyPoints> point_cloud;
//
//    ... Populate vector with points
//
//    hpm::ConvexHull2d<MyPoints> my_hull(point_cloud);
//    if (!my_hull.isValid()) {
//
//      ... error handling
//
//    }
//    // Get a copy of points to manipulate
//    std::vector<MyPoints> hull_points = my_hull.getHull();
//
//    // Or iterate points without storing a copy
//    for (auto i : my_hull.getHull()) {
//      std::cout << i.x() << ' ' << i.y() << '\n';
//    }
//

#pragma once
#include <vector>
#include <utility>
#include <algorithm>
#include <tuple>
#include <limits>

namespace hpm {
template <typename T>
class ConvexHull2d
{
public:
  ConvexHull2d() = default;
  explicit ConvexHull2d(std::vector<T> in);
  const std::vector<T>& getHull() const { return m_hull; };
  bool isValid() const { return m_is_valid; }
  bool isEmpty() const { return m_hull.empty(); }
  const std::pair<double, double>& getCentroid() { return m_centroid; };
  double getArea() { return m_area; };
  std::vector<std::pair<double, double>> rotateVertices(double angle) const; // Rotates around origin (0,0)
  std::vector<std::pair<double, double>> rotateVertices(double angle, const T& center) const;
  std::vector<std::pair<double, double>> rotateVertices(double angle, const std::pair<double, double>& center) const;
private:
  std::vector<T> m_hull;
  void computeCentroidAndArea();
  std::pair<double, double> m_centroid;
  double m_area;
  bool m_is_valid;
};

template <typename T>
ConvexHull2d<T>::ConvexHull2d(std::vector<T> in) :
  m_area(0.0),
  m_is_valid(false)
{
  const size_t kValidHullSize = 3;
  std::sort(in.begin(), in.end(),
    [](const T& lhs, const T& rhs) { return std::tie(lhs.x(), lhs.y()) < std::tie(rhs.x(), rhs.y()); });
  auto unique_end = std::unique(in.begin(), in.end(),
    [](const T& lhs, const T& rhs) { return std::tie(lhs.x(), lhs.y()) == std::tie(rhs.x(), rhs.y()); });
  in.erase(unique_end, in.end());
  size_t num_points = in.size();
  if (num_points < kValidHullSize) {
    m_hull = in;
    computeCentroidAndArea();
    return;
  }
  // Determines CW or CCW turn. Casting to double to allow greater range of coordinates for
  // integers without overflowing.
  auto cross = [](const T& lhs_start, const T& lhs_end, const T& rhs_start, const T& rhs_end) {
    double lhs_x = static_cast<double>(lhs_end.x()) - static_cast<double>(lhs_start.x());
    double lhs_y = static_cast<double>(lhs_end.y()) - static_cast<double>(lhs_start.y());
    double rhs_x = static_cast<double>(rhs_end.x()) - static_cast<double>(rhs_start.x());
    double rhs_y = static_cast<double>(rhs_end.y()) - static_cast<double>(rhs_start.y());
    return lhs_x * rhs_y - lhs_y * rhs_x;
  };

  // Andrew's Monotone Chain algorithm
  // Build lower hull. 
  m_hull.resize(2 * num_points);
  size_t j = 0;
  for (size_t i = 0; i < num_points; ++i) {
    while (j >= 2 &&
      cross(m_hull[j - 2], m_hull[j - 1], m_hull[j - 2], in[i]) <= 0) {
      --j;
    }
    m_hull[j++] = in[i];
  }

  // Build upper hull.
  size_t k = j + 1;   // Start index of upper hull
  for (size_t i = num_points - 1; i-- > 0;) {
    while (j >= k &&
      cross(m_hull[j - 2], m_hull[j - 1], m_hull[j - 2], in[i]) <= 0) {
      --j;
    }
    m_hull[j++] = in[i];
  }
  m_hull.resize(j - 1);   // Last point is the same as the first point, so remove it.

  computeCentroidAndArea();
  auto is_finite = [](const double& x) {
    return (x <= std::numeric_limits<double>::max() && x >= -std::numeric_limits<double>::max());
  };
  if (!is_finite(m_area) || !is_finite(m_centroid.first) || !is_finite(m_centroid.second)) {
    m_area = 0.0;
    m_hull.clear();
    return;
  }
  m_is_valid = true;
}

template <typename T>
void ConvexHull2d<T>::computeCentroidAndArea()
{
  const size_t num_vertices = m_hull.size();
  switch (num_vertices) {
  case 0:
    return;
  case 1:
    m_centroid.first = m_hull.front().x();
    m_centroid.second = m_hull.front().y();
    return;
  case 2:
    m_centroid.first = (m_hull.front().x() + m_hull.back().x()) / 2.0;
    m_centroid.second = (m_hull.front().y() + m_hull.back().y()) / 2.0;
    return;
  default:
    break;    // Continue if hull has at least 3 vertices.
  }

  double current_x = 0.0;
  double current_y = 0.0;
  double previous_x = 0.0;
  double previous_y = 0.0;
  double signed_area = 0.0;
  double area = 0.0;
  size_t j = num_vertices - 1;
  for (size_t i = 0; i < num_vertices; ++i) {
    current_x = m_hull[i].x();
    current_y = m_hull[i].y();
    previous_x = m_hull[j].x();
    previous_y = m_hull[j].y();
    area = current_x * previous_y - previous_x * current_y;
    j = i;
    signed_area += area;
    m_centroid.first += (current_x + previous_x) * area;
    m_centroid.second += (current_y + previous_y) * area;
  }
  signed_area *= 0.5;
  m_centroid.first /= (6.0 * signed_area);
  m_centroid.second /= (6.0 * signed_area);
  m_area = std::abs(signed_area);
  return;
}


template <typename T>
std::vector<std::pair<double, double>> ConvexHull2d<T>::rotateVertices(double angle) const
{
  return rotateVertices(angle, std::make_pair(0.0, 0.0));
}


template <typename T>
std::vector<std::pair<double, double>> ConvexHull2d<T>::rotateVertices(double angle, const T& center) const
{
  return rotateVertices(angle, std::make_pair(static_cast<double>(center.x()), static_cast<double>(center.y())));
}


template <typename T>
std::vector<std::pair<double, double>> ConvexHull2d<T>::rotateVertices(double angle, const std::pair<double, double>& center) const
{
  double x, y;
  double cosine = std::cos(angle);
  double sine = std::sin(angle);
  std::vector<std::pair<double, double>> out;
  for (auto i : m_hull) {
    x = (i.x() - center.first) * cosine - (i.y() - center.second) * sine;
    y = (i.x() - center.first) * sine + (i.y() - center.second) * cosine;
    out.push_back(std::make_pair(x + center.first, y + center.second));
  }
  return out;
}

} // namespace ConvexHull
