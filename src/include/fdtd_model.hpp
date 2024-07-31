#ifndef __XFDTD_MODEL_FDTD_MODEL_HPP__
#define __XFDTD_MODEL_FDTD_MODEL_HPP__

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <string_view>
#include <xtensor.hpp>
#include <xtensor/xtensor_forward.hpp>

#include "geometry.hpp"
#include <xfdtd/common/index_task.h>
#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>

class FDTDModel {
public:
  template <typename T> using MyArray = xt::xarray<T>;
  using LengthType = xfdtd::Real;
  using SizeType = std::size_t;
  using InsidePointRecord = MyArray<std::multiset<SizeType>>;
  using Vector = xfdtd::Vector;

public:
  struct Info {
    SizeType num_node;
    SizeType num_element;
  };

  struct NodesSet {
    MyArray<LengthType> x;
    MyArray<LengthType> y;
    MyArray<LengthType> z;
  };

  struct TriangleElementSet {
    MyArray<SizeType> n0;
    MyArray<SizeType> n1;
    MyArray<SizeType> n2;
  };

  struct Domain {
    LengthType min_x;
    LengthType min_y;
    LengthType min_z;
    LengthType max_x;
    LengthType max_y;
    LengthType max_z;

    auto size() const {
      return Vector{max_x - min_x, max_y - min_y, max_z - min_z};
    }

    auto start() const { return Vector{min_x, min_y, min_z}; }

    auto end() const { return Vector{max_x, max_y, max_z}; }
  };

  using Box = xfdtd::IndexTask;

public:
  FDTDModel(LengthType dx, LengthType dy, LengthType dz)
      : _dx{dx}, _dy{dy}, _dz{dz} {}

  auto dx() const { return _dx; }

  auto dy() const { return _dy; }

  auto dz() const { return _dz; }

  auto info() const { return _info; }

  auto nodes() const { return _nodes; }

  auto elements() const { return _elements; }

  auto domain() const { return _domain; }

  auto box() const { return _box; }

  auto init(std::string_view data_dir) -> void {
    _info = readInfo((std::filesystem::path{data_dir} / "info.txt").string());
    _nodes = readNodes((std::filesystem::path{data_dir} / "nodes.txt").string(),
                       _info.num_node);
    _elements = readElements(
        (std::filesystem::path{data_dir} / "elements.txt").string(),
        _info.num_element);

    _domain = generateDomainSize(_nodes);
    _box = generateBox(_domain, _dx, _dy, _dz);

    auto total_num_node_x = _box.xRange().size() + 1;
    auto total_num_node_y = _box.yRange().size() + 1;
    auto total_num_node_z = _box.zRange().size() + 1;

    _point_inside_xy =
        InsidePointRecord::from_shape({total_num_node_x, total_num_node_y});
    _point_inside_yz =
        InsidePointRecord::from_shape({total_num_node_y, total_num_node_z});
    _point_inside_zx =
        InsidePointRecord::from_shape({total_num_node_z, total_num_node_x});
    for (auto &&p : _point_inside_xy) {
      p.clear();
    }

    for (auto &&p : _point_inside_yz) {
      p.clear();
    }

    for (auto &&p : _point_inside_zx) {
      p.clear();
    }

    _grid_node_object = MyArray<int>::from_shape(
        {total_num_node_x, total_num_node_y, total_num_node_z});

    for (auto &&p : _grid_node_object) {
      p = 0;
    }
  }

  auto read(std::string_view data_dir) -> void {
    init(data_dir);

    const auto dx = this->dx();
    const auto dy = this->dy();
    const auto dz = this->dz();
    const auto domain = this->domain();
    auto &&grid_node = _grid_node_object;

    for (auto e{0}; e < _info.num_element; ++e) {
      auto p0 = Vector{_nodes.x(_elements.n0(e)), _nodes.y(_elements.n0(e)),
                       _nodes.z(_elements.n0(e))};
      auto p1 = Vector{_nodes.x(_elements.n1(e)), _nodes.y(_elements.n1(e)),
                       _nodes.z(_elements.n1(e))};
      auto p2 = Vector{_nodes.x(_elements.n2(e)), _nodes.y(_elements.n2(e)),
                       _nodes.z(_elements.n2(e))};
      gridZLineCrossElement(domain, p0, p1, p2, dx, dy, dz, grid_node,
                            _point_inside_xy);
      gridYLineCrossElement(domain, p0, p1, p2, dx, dy, dz, grid_node,
                            _point_inside_zx);
      gridXLineCrossElement(domain, p0, p1, p2, dx, dy, dz, grid_node,
                            _point_inside_yz);
    }

    const auto box = this->box();
    fillZAxis(box, _point_inside_xy, grid_node);
    fillYAxis(box, _point_inside_zx, grid_node);
    fillXAxis(box, _point_inside_yz, grid_node);

    // clear the memory
    for (auto &&p : _point_inside_xy) {
      p.clear();
    }

    for (auto &&p : _point_inside_yz) {
      p.clear();
    }

    for (auto &&p : _point_inside_zx) {
      p.clear();
    }
  }

  auto saveToMsh(std::string_view filename) -> void {
    saveToMsh(filename, _grid_node_object, _box, _dx, _dy, _dz, 1);
  }

  auto isInside(LengthType x, LengthType y, LengthType z,
                LengthType eps) const -> bool {
    const auto &domain_region = this->_domain;

    if ((x < domain_region.min_x - eps) || (domain_region.max_x) + eps < x) {
      return false;
    }
    if ((y < domain_region.min_y - eps) || (domain_region.max_y + eps < y)) {
      return false;
    }
    if ((z < domain_region.min_z - eps) || (domain_region.max_z + eps < z)) {
      return false;
    }

    const auto start = domain_region.start();
    auto i = static_cast<SizeType>(std::round((x - start.x()) / _dx));
    auto j = static_cast<SizeType>(std::round((y - start.y()) / _dy));
    auto k = static_cast<SizeType>(std::round((z - start.z()) / _dz));

    const auto &ob = this->_grid_node_object;

    if (ob(i, j, k) == 0 || ob(i + 1, j, k) == 0 || ob(i, j + 1, k) == 0 ||
        ob(i + 1, j + 1, k) == 0 || ob(i, j, k + 1) == 0 ||
        ob(i + 1, j, k + 1) == 0 || ob(i, j + 1, k + 1) == 0 ||
        ob(i + 1, j + 1, k + 1) == 0) {
      return false;
    }

    return true;
  }

  auto toString() const {
    std::stringstream ss;

    ss << "Resolution: " << _dx << " " << _dy << " " << _dz << "\n";
    ss << "Domain: " << "Start: " << _domain.start().toString()
       << " End: " << _domain.end().toString()
       << " Size: " << _domain.size().toString() << "\n";
    ss << "Box: " << _box.toString();
    return ss.str();
  }

private:
  LengthType _dx;
  LengthType _dy;
  LengthType _dz;

  Info _info;
  NodesSet _nodes;
  TriangleElementSet _elements;
  Domain _domain;
  Box _box;

  InsidePointRecord _point_inside_xy;
  InsidePointRecord _point_inside_yz;
  InsidePointRecord _point_inside_zx;

  MyArray<int> _grid_node_object;

  auto readInfo(std::string_view filename) -> Info {
    std::ifstream info{filename.data(), std::ios::in};
    if (!info.is_open()) {
      std::stringstream ss;
      ss << "readInfo: Failed to open file " << filename;
      throw std::runtime_error(ss.str());
    }

    SizeType num_node = 0;
    SizeType num_element = 0;
    info >> num_node >> num_element;
    info.close();

    return {num_node, num_element};
  }

  auto readNodes(std::string_view filename, SizeType num_node) -> NodesSet {
    std::ifstream nodes{filename.data(), std::ios::in};
    if (!nodes.is_open()) {
      std::stringstream ss;
      ss << "readNodes: Failed to open file " << filename;
      throw std::runtime_error(ss.str());
    }

    auto x = MyArray<LengthType>::from_shape({num_node});
    auto y = MyArray<LengthType>::from_shape({num_node});
    auto z = MyArray<LengthType>::from_shape({num_node});

    for (SizeType i = 0; i < num_node; ++i) {
      nodes >> x(i) >> y(i) >> z(i);
    }
    nodes.close();

    return {x, y, z};
  }

  auto readElements(std::string_view filename,
                    SizeType num_element) -> TriangleElementSet {
    std::ifstream elements{filename.data(), std::ios::in};
    if (!elements.is_open()) {
      std::stringstream ss;
      ss << "readElements: Failed to open file " << filename;
      throw std::runtime_error(ss.str());
    }
    auto n0 = MyArray<SizeType>::from_shape({num_element});
    auto n1 = MyArray<SizeType>::from_shape({num_element});
    auto n2 = MyArray<SizeType>::from_shape({num_element});

    for (SizeType i = 0; i < num_element && !elements.eof(); ++i) {
      elements >> n0(i) >> n1(i) >> n2(i);
    }
    elements.close();

    if (n0.size() != num_element || n1.size() != num_element ||
        n2.size() != num_element) {
      std::stringstream ss;
      ss << "readElements: The number of elements is not correct";
      throw std::runtime_error(ss.str());
    }

    return {n0, n1, n2};
  }

  auto generateDomainSize(const NodesSet &nodes) const -> Domain {
    auto domain = Domain{};
    domain.min_x = std::ranges::min(nodes.x);
    domain.min_y = std::ranges::min(nodes.y);
    domain.min_z = std::ranges::min(nodes.z);
    domain.max_x = std::ranges::max(nodes.x);
    domain.max_y = std::ranges::max(nodes.y);
    domain.max_z = std::ranges::max(nodes.z);
    return domain;
  }

  auto generateBox(Domain &domain, LengthType dx, LengthType dy,
                   LengthType dz) -> Box {
    auto size_x = domain.max_x - domain.min_x;
    auto size_y = domain.max_y - domain.min_y;
    auto size_z = domain.max_z - domain.min_z;

    auto num_x = static_cast<SizeType>(std::round(size_x / dx));
    auto num_y = static_cast<SizeType>(std::round(size_y / dy));
    auto num_z = static_cast<SizeType>(std::round(size_z / dz));

    return xfdtd::makeIndexTask(xfdtd::makeIndexRange(0, num_x),
                                xfdtd::makeIndexRange(0, num_y),
                                xfdtd::makeIndexRange(0, num_z));
  }

  auto generateBox(const Domain &domain, LengthType dx, LengthType dy,
                   LengthType dz, const Vector &p0, const Vector &p1,
                   const Vector &p2) -> Box {
    auto min_x = std::min({p0.x(), p1.x(), p2.x()});
    auto min_y = std::min({p0.y(), p1.y(), p2.y()});
    auto min_z = std::min({p0.z(), p1.z(), p2.z()});
    auto max_x = std::max({p0.x(), p1.x(), p2.x()});
    auto max_y = std::max({p0.y(), p1.y(), p2.y()});
    auto max_z = std::max({p0.z(), p1.z(), p2.z()});

    if (min_x < domain.min_x || min_y < domain.min_y || min_z < domain.min_z ||
        max_x > domain.max_x || max_y > domain.max_y || max_z > domain.max_z) {
      return {};
    }

    auto is = static_cast<SizeType>(std::round((min_x - domain.min_x) / dx));
    auto js = static_cast<SizeType>(std::round((min_y - domain.min_y) / dy));
    auto ks = static_cast<SizeType>(std::round((min_z - domain.min_z) / dz));
    auto ie = static_cast<SizeType>(std::round((max_x - domain.min_x) / dx));
    auto je = static_cast<SizeType>(std::round((max_y - domain.min_y) / dy));
    auto ke = static_cast<SizeType>(std::round((max_z - domain.min_z) / dz));

    return xfdtd::makeIndexTask(xfdtd::makeIndexRange(is, ie),
                                xfdtd::makeIndexRange(js, je),
                                xfdtd::makeIndexRange(ks, ke));
  };

  auto gridXLineCrossElement(const Domain &domain, const Vector &p0,
                             const Vector &p1, const Vector &p2, LengthType dx,
                             LengthType dy, LengthType dz, MyArray<int> &ob,
                             InsidePointRecord &point_inside) -> void {
    const auto start = domain.start();

    auto box = generateBox(_domain, dx, dy, dz, p0, p1, p2);

    auto plane = Plane<Vector, LengthType>{p0, p1, p2};
    auto [a, b, c, d] = plane.equation();

    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    for (SizeType j{js}; j <= je; ++j) {
      for (SizeType k{ks}; k <= ke; ++k) {
        auto y = j * dy + start.y();
        auto z = k * dz + start.z();
        auto x = (-b * y - c * z - d) / a;
        auto point = Vector{x, y, z};
        if (pointInTriangle(p0, p1, p2, point) &&
            !plane.parallel(Vector{1, 0, 0})) {
          auto i = static_cast<SizeType>(std::round((x - start.x()) / dx));
          ob.at(i, j, k) = 1;
          point_inside(j, k).insert({i});
        }
      }
    }
  }

  auto gridYLineCrossElement(const Domain &domain, const Vector &p0,
                             const Vector &p1, const Vector &p2, LengthType dx,
                             LengthType dy, LengthType dz, MyArray<int> &ob,
                             InsidePointRecord &point_inside) -> void {
    const auto start = domain.start();
    auto box = generateBox(_domain, dx, dy, dz, p0, p1, p2);

    auto plane = Plane<Vector, LengthType>{p0, p1, p2};
    auto [a, b, c, d] = plane.equation();

    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    for (SizeType k{ks}; k <= ke; ++k) {
      for (SizeType i{is}; i <= ie; ++i) {
        auto z = k * dz + start.z();
        auto x = i * dx + start.x();
        auto y = (-a * x - c * z - d) / b;
        auto point = Vector{x, y, z};
        if (pointInTriangle(p0, p1, p2, point) &&
            !plane.parallel(Vector{0, 1, 0})) {
          auto j = static_cast<SizeType>(std::round((y - start.y()) / dy));
          ob.at(i, j, k) = 1;
          point_inside(k, i).insert({j});
        }
      }
    }
  }

  auto gridZLineCrossElement(const Domain &domain, const Vector &p0,
                             const Vector &p1, const Vector &p2, LengthType dx,
                             LengthType dy, LengthType dz,
                             MyArray<int> &grid_node,
                             InsidePointRecord &point_inside) -> void {
    const auto start = domain.start();
    auto box = generateBox(_domain, dx, dy, dz, p0, p1, p2);

    auto plane = Plane<Vector, LengthType>{p0, p1, p2};
    auto [a, b, c, d] = plane.equation();

    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    for (SizeType i{is}; i <= ie; ++i) {
      for (SizeType j{js}; j <= je; ++j) {
        auto x = i * dx + start.x();
        auto y = j * dy + start.y();
        auto z = -(a * x + b * y + d) / c;
        auto point = Vector{x, y, z};
        if (pointInTriangle(p0, p1, p2, point) &&
            !plane.parallel(Vector{0, 0, 1})) {
          auto k = static_cast<SizeType>(std::round((z - start.z()) / dz));
          grid_node.at(i, j, k) = 1;
          point_inside(i, j).insert({k});
        }
      }
    }
  }

  auto fillXAxis(const Box &box, const InsidePointRecord &point_inside,
                 MyArray<int> &ob) -> void {
    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    for (SizeType j{js}; j <= je; ++j) {
      for (SizeType k{ks}; k <= ke; ++k) {
        if (point_inside(j, k).size() <= 1) {
          continue;
        }

        const auto &points = point_inside(j, k);
        auto p = points.begin();
        for (; p != points.end();) {
          auto q = p;
          std::advance(q, 1);
          if (q == points.end()) {
            break;
          }

          for (SizeType i{*p}; i < *q + 1; ++i) {
            ob(i, j, k) = 1;
          }

          std::advance(p, 2);
        }
      }
    }
  }

  auto fillYAxis(const Box &box, const InsidePointRecord &point_inside,
                 MyArray<int> &ob) -> void {
    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    for (SizeType k{ks}; k <= ke; ++k) {
      for (SizeType i{is}; i <= ie; ++i) {
        if (point_inside(k, i).size() <= 1) {
          continue;
        }

        auto &points = point_inside(k, i);
        auto p = points.begin();
        for (; p != points.end();) {
          auto q = p;
          std::advance(q, 1);
          if (q == points.end()) {
            break;
          }

          for (SizeType j{*p}; j < *q + 1; ++j) {
            ob(i, j, k) = 1;
          }

          std::advance(p, 2);
        }
      }
    }
  }

  auto fillZAxis(const Box &box, const InsidePointRecord &point_inside,
                 MyArray<int> &ob) -> void {
    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    for (SizeType i{is}; i <= ie; ++i) {
      for (SizeType j{js}; j <= je; ++j) {
        if (point_inside(i, j).size() <= 1) {
          continue;
        }

        auto &points = point_inside(i, j);

        auto p = points.begin();
        for (; p != points.end();) {
          auto q = p;
          std::advance(q, 1);
          if (q == points.end()) {
            break;
          }

          for (SizeType k{*p}; k < *q + 1; ++k) {
            ob(i, j, k) = 1;
          }

          std::advance(p, 2);
        }
      }
    }
  }

  auto saveToMsh(std::string_view path, const MyArray<int> &grid_node_object,
                 const Box &box, LengthType dx, LengthType dy, LengthType dz,
                 int object_index) const -> void {
    if (!box.valid()) {
      std::stringstream ss;
      ss << "saveToMsh: Invalid box";
      ss << box.toString();
      throw std::runtime_error(ss.str());
    }

    std::ofstream msh(path.data(), std::ios::out | std::ios::trunc);
    if (!msh.is_open()) {
      std::stringstream ss;
      ss << "saveToMsh: Failed to open file " << path;
      throw std::runtime_error(ss.str());
    }

    const auto is = box.xRange().start();
    const auto js = box.yRange().start();
    const auto ks = box.zRange().start();
    const auto ie = box.xRange().end();
    const auto je = box.yRange().end();
    const auto ke = box.zRange().end();

    auto nodes = MyArray<SizeType>::from_shape(grid_node_object.shape());

    msh << "$MeshFormat\n";
    msh << "2.2 0 8\n";
    msh << "$EndMeshFormat\n";

    {
      msh << "$Nodes\n";

      SizeType num_fdtd_grid = std::accumulate(
          grid_node_object.begin(), grid_node_object.end(), 0,
          [&](int acc, int v) { return acc + (v == object_index ? 1 : 0); });

      msh << num_fdtd_grid << "\n";

      SizeType counter = 0;
      for (auto i{is}; i <= ie; ++i) {
        for (auto j{js}; j <= je; ++j) {
          for (auto k{ks}; k <= ke; ++k) {
            if (grid_node_object(i, j, k) != object_index) {
              nodes(i, j, k) = 0;
              continue;
            }

            counter++;
            msh << counter << " " << domain().start().x() + i * dx << " "
                << domain().start().y() + j * dy << " "
                << domain().start().z() + k * dz << "\n";
            nodes(i, j, k) = counter;
          }
        }
      }

      msh << "$EndNodes\n";
    }

    {
      msh << "$Elements\n";

      SizeType num_fdtd_element = 0;

      for (auto i{is}; i < ie; ++i) {
        for (auto j{js}; j < je; ++j) {
          for (auto k{ks}; k < ke; ++k) {
            if (!gridContainObject(grid_node_object, i, j, k, object_index)) {
              continue;
            }
            ++num_fdtd_element;
          }
        }
      }

      msh << num_fdtd_element << "\n";

      SizeType counter = 0;
      for (auto i{is}; i < ie; ++i) {
        for (auto j{js}; j < je; ++j) {
          for (auto k{ks}; k < ke; ++k) {
            if (!gridContainObject(grid_node_object, i, j, k, object_index)) {
              continue;
            }

            counter++;
            msh << counter << " " << 5 << " " << 2 << " "
                << grid_node_object(i, j, k) << " "
                << grid_node_object(i + 1, j, k) << " " << nodes(i, j, k) << " "
                << nodes(i, j, k) << " " << nodes(i + 1, j + 1, k) << " "
                << nodes(i, j + 1, k) << " " << nodes(i, j, k + 1) << " "
                << nodes(i + 1, j, k + 1) << " " << nodes(i + 1, j + 1, k + 1)
                << " " << nodes(i, j + 1, k + 1) << "\n";
          }
        }
      }

      msh << "$EndElements";
    }

    msh.close();
  }

  auto gridContainObject(const MyArray<int> &grid_node_object, SizeType i,
                         SizeType j, SizeType k,
                         int object_index) const -> bool {
    const auto point_0 = grid_node_object(i, j, k);
    const auto point_1 = grid_node_object(i + 1, j, k);
    const auto point_2 = grid_node_object(i, j + 1, k);
    const auto point_3 = grid_node_object(i + 1, j + 1, k);
    const auto point_4 = grid_node_object(i, j, k + 1);
    const auto point_5 = grid_node_object(i + 1, j, k + 1);
    const auto point_6 = grid_node_object(i, j + 1, k + 1);
    const auto point_7 = grid_node_object(i + 1, j + 1, k + 1);

    return point_0 == object_index && point_1 == object_index &&
           point_2 == object_index && point_3 == object_index &&
           point_4 == object_index && point_5 == object_index &&
           point_6 == object_index && point_7 == object_index;
  }
};

#endif // __XFDTD_MODEL_FDTD_MODEL_HPP__
