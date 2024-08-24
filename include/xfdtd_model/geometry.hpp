#ifndef __XFDTD_MODEL_SUPPORT_GEOMETRY_HPP__
#define __XFDTD_MODEL_SUPPORT_GEOMETRY_HPP__

#include <xfdtd/coordinate_system/coordinate_system.h>

#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <sstream>

template <typename Vector = xfdtd::Vector, typename T = xfdtd::Real>
class Plane {
 public:
  struct PlaneEquation {
    T _a, _b, _c, _d;
  };

  static auto calculateEquation(const Vector &p0, const Vector &p1,
                                const Vector &p2) -> PlaneEquation;

 public:
  Plane(Vector p0, Vector p1, Vector p2);

  Plane(const Vector &normal, const Vector &p);

  auto equation() const -> PlaneEquation;

  auto normal() const -> Vector;

  auto parallel(const Vector &v) const -> bool;

  auto pointInPlane(const Vector &p) const -> bool;

  auto p0() const -> Vector { return _p0; }

  auto p1() const -> Vector { return _p1; }

  auto p2() const -> Vector { return _p2; }

 private:
  Vector _p0{}, _p1{}, _p2{};
  PlaneEquation _eq{};
};

template <typename Vector, typename T>
inline auto Plane<Vector, T>::calculateEquation(
    const Vector &p0, const Vector &p1, const Vector &p2) -> PlaneEquation {
  auto &&v = p1 - p0;
  auto &&u = p2 - p0;
  auto tem = v - u;
  // auto l = tem.normL2();
  auto &&n = v.cross(u);
  n = n / n.normL2();

  auto a = n.x();
  auto b = n.y();
  auto c = n.z();
  auto d = -n.dot(p0);
  return {a, b, c, d};
}

template <typename Vector, typename T>
inline Plane<Vector, T>::Plane(Vector p0, Vector p1, Vector p2)
    : _p0{p0}, _p1{p1}, _p2{p2}, _eq{calculateEquation(p0, p1, p2)} {}

template <typename Vector, typename T>
inline Plane<Vector, T>::Plane(const Vector &normal, const Vector &p) {
  // check if the normal is zero
  if (normal.normL2() < Vector::zero_tolerance) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6);
    ss << "normal: " << normal << "\np: " << p << "\n";
    throw std::invalid_argument("The normal is zero: \n" + ss.str());
    // std::cerr << "The normal is zero: \n" + ss.str();
  }

  // check if the point is in the normal
  if (normal.orthogonal(p)) {
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6);
    ss << "normal: " << normal << "\np: " << p << "\n";
    throw std::invalid_argument("The point is in the normal: \n" + ss.str());
    // std::cerr << "The point is in the normal: \n" + ss.str();
  }

  auto d = -normal.dot(p);
  auto a = normal.x();
  auto b = normal.y();
  auto c = normal.z();
  _p0 = p;
  _p1 = p + Vector(1, 0, 0);
  _p2 = p + Vector(0, 1, 0);
  _eq = {a, b, c, d};
}

template <typename Vector, typename T>
inline auto Plane<Vector, T>::equation() const -> PlaneEquation {
  return calculateEquation(_p0, _p1, _p2);
}

template <typename Vector, typename T>
inline auto Plane<Vector, T>::normal() const -> Vector {
  return Vector(_eq._a, _eq._b, _eq._c);
}

template <typename Vector, typename T>
inline auto Plane<Vector, T>::parallel(const Vector &v) const -> bool {
  return normal().orthogonal(v);
}

template <typename Vector, bool contain_edge = true>
inline auto pointInTriangle(const Vector &a, const Vector &b, const Vector &c,
                            const Vector &p) -> bool {
  auto &&v0 = c - a;
  auto &&v1 = b - a;
  auto &&v2 = p - a;

  auto dot_00 = v0.dot(v0);
  auto dot_01 = v0.dot(v1);
  auto dot_02 = v0.dot(v2);
  auto dot_11 = v1.dot(v1);
  auto dot_12 = v1.dot(v2);

  auto inverse = 1 / (dot_00 * dot_11 - dot_01 * dot_01);

  auto u = (dot_11 * dot_02 - dot_01 * dot_12) * inverse;
  auto v = (dot_00 * dot_12 - dot_01 * dot_02) * inverse;

  if constexpr (contain_edge) {
    if (u < 0 || u > 1) {
      return false;
    }

    if (v < 0 || v > 1) {
      return false;
    }

    return u + v <= 1;
  } else {
    if (u <= 0 || u >= 1) {
      // u is (0, 1). If u is 0 or 1, the point is on the edge of the triangle.
      return false;
    }

    if (v <= 0 || v >= 1) {
      return false;
    }

    return u + v < 1;
  }
}

#endif  // __XFDTD_MODEL_SUPPORT_GEOMETRY_HPP__
