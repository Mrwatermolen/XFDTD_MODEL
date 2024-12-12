#ifndef __XFDTD_MODEL_GRID_MODEL_HPP__
#define __XFDTD_MODEL_GRID_MODEL_HPP__

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>
#include <xfdtd/util/transform/abc_xyz.h>

#include <iterator>
#include <set>
#include <sstream>
#include <utility>
#include <xfdtd_model/geometry.hpp>
#include <xfdtd_model/triangular_model_info.hpp>
#include <xtensor.hpp>

namespace xfdtd::model {

class GridModel {
 public:
  using IntersectionRecord = Array2D<std::set<xfdtd::Index>>;

  auto& triangularModelInfo() const { return _triangular_model_info; }

  auto& triangularModelInfo() { return _triangular_model_info; }

  auto read(std::istream& info_is, std::istream& vertices_is,
            std::istream& faces_is) -> void {
    _triangular_model_info =
        makeTriangularModelInfo(info_is, vertices_is, faces_is);
  }

  auto buildUniform(Real dx, Real dy, Real dz) -> void {
    auto&& triangular_model_info = _triangular_model_info;
    Array1D<Real> x_node_positions = xt::arange<Real>(
        triangular_model_info.minX(), triangular_model_info.maxX() + dx, dx);
    Array1D<Real> y_node_positions = xt::arange<Real>(
        triangular_model_info.minY(), triangular_model_info.maxY() + dy, dy);
    Array1D<Real> z_node_positions = xt::arange<Real>(
        triangular_model_info.minZ(), triangular_model_info.maxZ() + dz, dz);
    buildModel(x_node_positions, y_node_positions, z_node_positions);
  }

  auto buildModel(Array1D<Real> x_node_positions,
                  Array1D<Real> y_node_positions,
                  Array1D<Real> z_node_positions) -> void {
    _x_node_positions = std::move(x_node_positions);
    _y_node_positions = std::move(y_node_positions);
    _z_node_positions = std::move(z_node_positions);

    _node_record = buildNodeRecord(_triangular_model_info, _x_node_positions,
                                   _y_node_positions, _z_node_positions);

    _x_center_positions = (xt::view(_x_node_positions, xt::range(_, -1)) +
                           xt::view(_x_node_positions, xt::range(1, _))) *
                          0.5;
    _y_center_positions = (xt::view(_y_node_positions, xt::range(_, -1)) +
                           xt::view(_y_node_positions, xt::range(1, _))) *
                          0.5;
    _z_center_positions = (xt::view(_z_node_positions, xt::range(_, -1)) +
                           xt::view(_z_node_positions, xt::range(1, _))) *
                          0.5;
  }

  auto writeToMsh(std::ostream& msh) -> void {
    writeToMsh(msh, _node_record, _x_node_positions, _y_node_positions,
               _z_node_positions);
  }

  auto isInside(Real x, Real y, Real z, Real eps) const -> bool {
    if (x + eps < _x_node_positions.front() ||
        _x_node_positions.back() < x - eps) {
      return false;
    }

    if (y + eps < _y_node_positions.front() ||
        _y_node_positions.back() < y - eps) {
      return false;
    }

    if (z + eps < _z_node_positions.front() ||
        _z_node_positions.back() < z - eps) {
      return false;
    }

    auto i = xt::argmin(xt::abs(_x_center_positions - x)).front();
    auto j = xt::argmin(xt::abs(_y_center_positions - y)).front();
    auto k = xt::argmin(xt::abs(_z_center_positions - z)).front();
    return isObjectGrid(_node_record, i, j, k);
  }

 private:
  TriangularModelInfo _triangular_model_info{};
  Array3D<int> _node_record;

  Array1D<Real> _x_node_positions, _y_node_positions, _z_node_positions;
  Array1D<Real> _x_center_positions, _y_center_positions, _z_center_positions;

  static auto buildNodeRecord(const TriangularModelInfo& triangular_model_info,
                              const Array1D<Real>& x_node_positions,
                              const Array1D<Real>& y_node_positions,
                              const Array1D<Real>& z_node_positions)
      -> Array3D<int> {
    auto nx = x_node_positions.size();
    auto ny = y_node_positions.size();
    auto nz = z_node_positions.size();
    if (nx == 0 || ny == 0 || nz == 0) {
      std::stringstream ss;
      ss << "GridModel::build: Invalid grid size: nx=" << nx << ", ny=" << ny
         << ", nz=" << nz;
      throw std::runtime_error(ss.str());
    }

    {
      auto n = std::numeric_limits<Index>::max() / nx;
      if (n < ny * nz) {
        std::stringstream ss;
        ss << "GridModel::build: Too many grids: nx=" << nx << ", ny=" << ny
           << ", nz=" << nz;
        throw std::runtime_error(ss.str());
      }
    }

    auto intersection_record_yz = IntersectionRecord({ny, nz});
    auto intersection_record_zx = IntersectionRecord({nz, nx});
    auto intersection_record_xy = IntersectionRecord({nx, ny});
    for (auto&& i : intersection_record_yz) {
      i.clear();
    }
    for (auto&& i : intersection_record_zx) {
      i.clear();
    }
    for (auto&& i : intersection_record_xy) {
      i.clear();
    }

    auto node_record = buildNodeRecord(
        triangular_model_info, x_node_positions, y_node_positions,
        z_node_positions, intersection_record_yz, intersection_record_zx,
        intersection_record_xy);

    fill<xfdtd::Axis::XYZ::X>(x_node_positions, y_node_positions,
                              z_node_positions, intersection_record_yz,
                              node_record);
    fill<xfdtd::Axis::XYZ::Y>(x_node_positions, y_node_positions,
                              z_node_positions, intersection_record_zx,
                              node_record);
    fill<xfdtd::Axis::XYZ::Z>(x_node_positions, y_node_positions,
                              z_node_positions, intersection_record_xy,
                              node_record);

    return node_record;
  }

  static auto makeTriangularModelInfo(
      std::istream& info_is, std::istream& vertices_is,
      std::istream& faces_is) -> TriangularModelInfo {
    TriangularModelInfo triangular_model_info;
    triangular_model_info.read(info_is, vertices_is, faces_is);
    return triangular_model_info;
  }

  static auto buildNodeRecord(const TriangularModelInfo& triangular_model_info,
                              const Array1D<Real>& x_node_positions,
                              const Array1D<Real>& y_node_positions,
                              const Array1D<Real>& z_node_positions,
                              IntersectionRecord& intersection_record_yz,
                              IntersectionRecord& intersection_record_zx,
                              IntersectionRecord& intersection_record_xy)
      -> Array3D<int> {
    const auto nx = x_node_positions.size();
    const auto ny = y_node_positions.size();
    const auto nz = z_node_positions.size();

    Array3D<int> node_record({nx, ny, nz});
    node_record.fill(0);

    for (const auto& element : triangular_model_info.elements()) {
      auto p0 = triangular_model_info.vertices()[element.n0()];
      auto p1 = triangular_model_info.vertices()[element.n1()];
      auto p2 = triangular_model_info.vertices()[element.n2()];

      detectIntersection<xfdtd::Axis::XYZ::X>(
          p0, p1, p2, x_node_positions, y_node_positions, z_node_positions,
          node_record, intersection_record_yz);
      detectIntersection<xfdtd::Axis::XYZ::Y>(
          p0, p1, p2, x_node_positions, y_node_positions, z_node_positions,
          node_record, intersection_record_zx);
      detectIntersection<xfdtd::Axis::XYZ::Z>(
          p0, p1, p2, x_node_positions, y_node_positions, z_node_positions,
          node_record, intersection_record_xy);
    }

    return node_record;
  };

  static auto writeToMsh(std::ostream& msh, const Array3D<int>& node_record,
                         const Array1D<Real>& x_node_positions,
                         const Array1D<Real>& y_node_positions,
                         const Array1D<Real>& z_node_positions) -> void {
    const Index is = 0;
    const Index ie = node_record.shape()[0] - 1;
    const Index js = 0;
    const Index je = node_record.shape()[1] - 1;
    const Index ks = 0;
    const Index ke = node_record.shape()[2] - 1;

    auto node_index = Array3D<Index>({ie - is + 1, je - js + 1, ke - ks + 1});
    node_index.fill(0);

    msh << "$MeshFormat\n";
    msh << "2.2 0 8\n";
    msh << "$EndMeshFormat\n";

    {
      msh << "$Nodes\n";
      auto num_grid_in_object = std::accumulate(
          node_record.begin(), node_record.end(), Index{0},
          [&](auto acc, auto v) { return acc + (v == 1 ? 1 : 0); });

      msh << num_grid_in_object << '\n';

      Index node_index_count = 1;
      for (Index i{is}; i <= ie; ++i) {
        for (Index j{js}; j <= je; ++j) {
          for (Index k{ks}; k <= ke; ++k) {
            if (node_record(i, j, k) != 1) {
              continue;
            }

            node_index.at(i, j, k) = node_index_count;
            msh << node_index_count << ' ' << x_node_positions[i] << ' '
                << y_node_positions[j] << ' ' << z_node_positions[k] << '\n';
            node_index_count++;
          }
        }
      }

      msh << "$EndNodes\n";
    }

    {
      msh << "$Elements\n";
      Index num_elements = 0;

      for (auto i{is}; i < ie; ++i) {
        for (auto j{js}; j < je; ++j) {
          for (auto k{ks}; k < ke; ++k) {
            if (!isObjectGrid(node_record, i, j, k)) {
              continue;
            }

            ++num_elements;
          }
        }
      }

      msh << num_elements << "\n";

      Index element_index = 1;
      for (auto i{is}; i < ie; ++i) {
        for (auto j{js}; j < je; ++j) {
          for (auto k{ks}; k < ke; ++k) {
            if (!isObjectGrid(node_record, i, j, k)) {
              continue;
            }

            msh << element_index << " " << 5 << " " << 2 << " " << 1 << " " << 1
                << " " << node_index(i, j, k) << " " << node_index(i, j, k)
                << " " << node_index(i + 1, j + 1, k) << " "
                << node_index(i, j + 1, k) << " " << node_index(i, j, k + 1)
                << " " << node_index(i + 1, j, k + 1) << " "
                << node_index(i + 1, j + 1, k + 1) << " "
                << node_index(i, j + 1, k + 1) << "\n";
            element_index++;
          }
        }
      }

      msh << "$EndElements";
    }
  }

  template <xfdtd::Axis::XYZ xyz>
  static auto detectIntersection(
      const Vector& p0, const Vector& p1, const Vector& p2,
      const Array1D<Real>& x_node_positions,
      const Array1D<Real>& y_node_positions,
      const Array1D<Real>& z_node_positions, Array3D<int>& node_record,
      IntersectionRecord& intersection_record) -> void {
    auto plane = Plane{p0, p1, p2};
    auto [a, b, c, d] = plane.equation();

    Index is = 0;
    Index ie = 0;
    Index js = 0;
    Index je = 0;
    Index ks = 0;
    Index ke = 0;

    {
      auto triangle_min_x = std::min({p0.x(), p1.x(), p2.x()});
      auto triangle_max_x = std::max({p0.x(), p1.x(), p2.x()});
      auto triangle_min_y = std::min({p0.y(), p1.y(), p2.y()});
      auto triangle_max_y = std::max({p0.y(), p1.y(), p2.y()});
      auto triangle_min_z = std::min({p0.z(), p1.z(), p2.z()});
      auto triangle_max_z = std::max({p0.z(), p1.z(), p2.z()});
      auto get_floor_index = [](const auto& positions,
                                const auto& triangle_min) {
        auto index = xt::argmin(xt::abs(positions - triangle_min)).front();
        if (triangle_min < positions[index] && index != 0) {
          index -= 1;
        }
        return index;
      };
      auto get_ceil_index = [](const auto& positions,
                               const auto& triangle_max) {
        auto index = xt::argmin(xt::abs(positions - triangle_max)).front();
        if (positions[index] < triangle_max && index < positions.size() - 1) {
          index += 1;
        }
        return index;
      };

      is = get_floor_index(x_node_positions, triangle_min_x);
      ie = std::min(get_ceil_index(x_node_positions, triangle_max_x) + 1,
                    x_node_positions.size());
      js = get_floor_index(y_node_positions, triangle_min_y);
      je = std::min(get_ceil_index(y_node_positions, triangle_max_y) + 1,
                    y_node_positions.size());
      ks = get_floor_index(z_node_positions, triangle_min_z);
      ke = std::min(get_ceil_index(z_node_positions, triangle_max_z) + 1,
                    z_node_positions.size());
    }

    auto set_node_intersections = [&node_record](auto i, auto j, auto k) {
      node_record.at(i, j, k) = 1;
    };
    auto record_intersection = [&intersection_record](auto a, auto b, auto c) {
      intersection_record.at(a, b).insert(c);
    };

    if constexpr (xyz == xfdtd::Axis::XYZ::X) {
      auto get_index_func = [&x_node_positions](const auto& point) {
        return xt::argmin(xt::abs(x_node_positions - point)).front();
      };

      auto normal = Vector{1, 0, 0};

      for (auto j{js}; j < je; ++j) {
        for (auto k{ks}; k < ke; ++k) {
          const auto y = y_node_positions[j];
          const auto z = z_node_positions[k];
          const auto x = -(b * y + c * z + d) / a;

          auto p = Vector{x, y, z};
          if (pointInTriangle(p0, p1, p2, p) && !plane.parallel(normal)) {
            // Intersection
            auto i = get_index_func(x);
            set_node_intersections(i, j, k);
            auto [a, b, c] = xfdtd::transform::xYZToABC<Index, xyz>(i, j, k);
            record_intersection(a, b, c);
          }
        }
      }
    }

    if constexpr (xyz == xfdtd::Axis::XYZ::Y) {
      auto get_index_func = [&y_node_positions](const auto& point) {
        return xt::argmin(xt::abs(y_node_positions - point)).front();
      };

      auto normal = Vector{0, 1, 0};

      for (auto i{is}; i < ie; ++i) {
        for (auto k{ks}; k < ke; ++k) {
          const auto x = x_node_positions[i];
          const auto z = z_node_positions[k];
          const auto y = -(a * x + c * z + d) / b;

          auto p = Vector{x, y, z};
          if (pointInTriangle(p0, p1, p2, p) && !plane.parallel(normal)) {
            // Intersection
            auto j = get_index_func(y);
            set_node_intersections(i, j, k);
            auto [a, b, c] = xfdtd::transform::xYZToABC<Index, xyz>(i, j, k);
            record_intersection(a, b, c);
          }
        }
      }
    }

    if constexpr (xyz == xfdtd::Axis::XYZ::Z) {
      auto get_index_func = [&z_node_positions](const auto& point) {
        return xt::argmin(xt::abs(z_node_positions - point)).front();
      };

      auto normal = Vector{0, 0, 1};

      for (auto i{is}; i < ie; ++i) {
        for (auto j{js}; j < je; ++j) {
          const auto x = x_node_positions[i];
          const auto y = y_node_positions[j];
          const auto z = -(a * x + b * y + d) / c;

          auto p = Vector{x, y, z};
          if (pointInTriangle(p0, p1, p2, p) && !plane.parallel(normal)) {
            // Intersection
            auto k = get_index_func(z);
            set_node_intersections(i, j, k);
            auto [a, b, c] = xfdtd::transform::xYZToABC<Index, xyz>(i, j, k);
            record_intersection(a, b, c);
          }
        }
      }
    }
  }

  template <xfdtd::Axis::XYZ xyz>
  static auto fill(const Array1D<Real>& x_node_positions,
                   const Array1D<Real>& y_node_positions,
                   const Array1D<Real>& z_node_positions,
                   const IntersectionRecord& intersection_record,
                   Array3D<int>& node_record) -> void {
    auto is = Index{0};
    auto ie = x_node_positions.size();
    auto js = Index{0};
    auto je = y_node_positions.size();
    auto ks = Index{0};
    auto ke = z_node_positions.size();

    {
      auto [a, b, c] = xfdtd::transform::xYZToABC<Index, xyz>(ie, je, ke);
      auto [i, j, k] = xfdtd::transform::aBCToXYZ<Index, xyz>(a, b, 1);
      ie = i;
      je = j;
      ke = k;
    }

    for (auto i{is}; i < ie; ++i) {
      for (auto j{js}; j < je; ++j) {
        for (auto k{ks}; k < ke; ++k) {
          auto [a, b, c] = xfdtd::transform::xYZToABC<Index, xyz>(i, j, k);
          const auto& points_index = intersection_record(a, b);
          if (points_index.empty()) {
            continue;
          }

          auto p = points_index.begin();
          while (p != points_index.end()) {
            auto q = p;
            std::advance(q, 1);
            if (q == points_index.end()) {
              break;
            }

            for (auto c0 = *p; c0 < *q + 1; ++c0) {
              auto [i0, j0, k0] =
                  xfdtd::transform::aBCToXYZ<Index, xyz>(a, b, c0);
              node_record.at(i0, j0, k0) = 1;
            }

            std::advance(p, 2);
          }
        }
      }
    }
  }

  static auto isObjectGrid(const Array3D<int>& node_record, Index i, Index j,
                           Index k) -> bool {
    const auto point_0 = node_record(i, j, k);
    const auto point_1 = node_record(i + 1, j, k);
    const auto point_2 = node_record(i, j + 1, k);
    const auto point_3 = node_record(i + 1, j + 1, k);
    const auto point_4 = node_record(i, j, k + 1);
    const auto point_5 = node_record(i + 1, j, k + 1);
    const auto point_6 = node_record(i, j + 1, k + 1);
    const auto point_7 = node_record(i + 1, j + 1, k + 1);

    return point_0 == 1 && point_1 == 1 && point_2 == 1 && point_3 == 1 &&
           point_4 == 1 && point_5 == 1 && point_6 == 1 && point_7 == 1;
  }
};

}  // namespace xfdtd::model

#endif  // __XFDTD_MODEL_GRID_MODEL_HPP__
