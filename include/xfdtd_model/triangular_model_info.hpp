#ifndef __XFDTD_MODEL_TRIANGULAR_MODEL_INFO_HPP__
#define __XFDTD_MODEL_TRIANGULAR_MODEL_INFO_HPP__

#include <xfdtd/common/type_define.h>
#include <xfdtd/coordinate_system/coordinate_system.h>

#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <sstream>
#include <vector>

namespace xfdtd::model {

class TriangularModelInfo {
 public:
  struct TriangularElement {
    Index _vertex_index[3]{};

    auto n0() const -> Index { return _vertex_index[0]; }

    auto n1() const -> Index { return _vertex_index[1]; }

    auto n2() const -> Index { return _vertex_index[2]; }
  };

  auto numVertices() const -> Index { return _num_vertices; }

  auto numElements() const -> Index { return _num_elements; }

  auto minX() const -> Real { return _min_x; }

  auto minY() const -> Real { return _min_y; }

  auto minZ() const -> Real { return _min_z; }

  auto sizeX() const -> Real { return _size_x; }

  auto sizeY() const -> Real { return _size_y; }

  auto sizeZ() const -> Real { return _size_z; }

  auto maxX() const -> Real { return _min_x + _size_x; }

  auto maxY() const -> Real { return _min_y + _size_y; }

  auto maxZ() const -> Real { return _min_z + _size_z; }

  auto& vertices() const { return _vertices; }

  auto& elements() const { return _elements; }

  auto read(std::istream& info_is, std::istream& vertices_is,
            std::istream& faces_is) -> void {
    if (!info_is) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Invalid input stream";
      throw std::runtime_error(ss.str());
    }

    _num_vertices = readNumVertices(info_is);

    if (!info_is && !info_is.eof()) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Failed to read number of vertices";
      throw std::runtime_error(ss.str());
    }

    if (_num_vertices == 0) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: No vertices found";
    }

    _num_elements = readNumFaces(info_is);

    if (!info_is && !info_is.eof()) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Failed to read number of faces";
      throw std::runtime_error(ss.str());
    }

    if (_num_elements == 0) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: No faces found";
      throw std::runtime_error(ss.str());
    }

    _vertices = readVertices(vertices_is);

    if (!vertices_is && !vertices_is.eof()) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Failed to read vertices";
      throw std::runtime_error(ss.str());
    }

    if (_vertices.size() != _num_vertices) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Number of vertices mismatch";
      ss << " (expected: " << _num_vertices << ", actual: " << _vertices.size()
         << ")";
      throw std::runtime_error(ss.str());
    }

    _elements = readElements(faces_is);

    if (!faces_is && !faces_is.eof()) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Failed to read faces";
      throw std::runtime_error(ss.str());
    }

    if (_elements.size() != _num_elements) {
      std::stringstream ss;
      ss << "TriangleObjectInfo::read: Number of faces mismatch";
      ss << " (expected: " << _num_elements << ", actual: " << _elements.size()
         << ")";
      throw std::runtime_error(ss.str());
    }

    _min_x = _min_y = _min_z = std::numeric_limits<Real>::max();
    Real max_x = std::numeric_limits<Real>::lowest();
    Real max_y = std::numeric_limits<Real>::lowest();
    Real max_z = std::numeric_limits<Real>::lowest();

    for (const auto& vertex : _vertices) {
      if (vertex.x() < _min_x) {
        _min_x = vertex.x();
      }

      if (vertex.y() < _min_y) {
        _min_y = vertex.y();
      }

      if (vertex.z() < _min_z) {
        _min_z = vertex.z();
      }

      if (vertex.x() > max_x) {
        max_x = vertex.x();
      }

      if (vertex.y() > max_y) {
        max_y = vertex.y();
      }

      if (vertex.z() > max_z) {
        max_z = vertex.z();
      }
    }

    _size_x = max_x - _min_x;
    _size_y = max_y - _min_y;
    _size_z = max_z - _min_z;

    if (_size_x < 0 || _size_y < 0 || _size_z < 0) {
      std::stringstream ss;
      ss << std::setprecision(16) << std::scientific;
      ss << "TriangleObjectInfo::read: Invalid size";
      ss << " (size_x: " << _size_x << ", size_y: " << _size_y
         << ", size_z: " << _size_z << ")";
      throw std::runtime_error(ss.str());
    }

    if (std::abs(_size_x) < std::numeric_limits<Real>::epsilon() ||
        std::abs(_size_y) < std::numeric_limits<Real>::epsilon() ||
        std::abs(_size_z) < std::numeric_limits<Real>::epsilon()) {
      std::stringstream ss;
      ss << std::setprecision(16) << std::scientific;
      ss << "TriangleObjectInfo::read: Invalid size";
      ss << " (size_x: " << _size_x << ", size_y: " << _size_y
         << ", size_z: " << _size_z << ")";
      throw std::runtime_error(ss.str());
    }
  }

 private:
  static auto readNumVertices(std::istream& is) -> Index {
    Index num_vertices{};
    is >> num_vertices;
    return num_vertices;
  }

  static auto readNumFaces(std::istream& is) -> Index {
    Index num_faces{};
    is >> num_faces;
    return num_faces;
  }

  static auto readVertices(std::istream& is) -> std::vector<Vector> {
    std::vector<Vector> vertices;
    while (is) {
      Real x{};
      Real y{};
      Real z{};
      is >> x >> y >> z;
      if (!is) {
        break;
      }

      vertices.emplace_back(x, y, z);
    }

    return vertices;
  }

  static auto readElements(std::istream& is) -> std::vector<TriangularElement> {
    std::vector<TriangularElement> elements;
    while (is) {
      TriangularElement element;
      is >> element._vertex_index[0] >> element._vertex_index[1] >>
          element._vertex_index[2];
      if (!is) {
        break;
      }

      elements.emplace_back(element);
    }

    return elements;
  }

 private:
  Index _num_vertices{};
  Index _num_elements{};

  Real _min_x{}, _min_y{}, _min_z{};
  Real _size_x{}, _size_y{}, _size_z{};
  std::vector<Vector> _vertices;
  std::vector<TriangularElement> _elements;
};

}  // namespace xfdtd::model

#endif  // __XFDTD_MODEL_TRIANGULAR_MODEL_INFO_H__
