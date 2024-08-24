#ifndef __XFDTD_MODEL_MODEL_SHAPE_HPP__
#define __XFDTD_MODEL_MODEL_SHAPE_HPP__

#include <xfdtd/common/type_define.h>
#include <xfdtd/shape/shape.h>

#include <istream>
#include <memory>
#include <utility>
#include <xfdtd_model/grid_model.hpp>

namespace xfdtd::model {

class ModeShapeBase : public Shape {
 public:
  ModeShapeBase() = default;

  ~ModeShapeBase() override = default;

  ModeShapeBase(const ModeShapeBase&) = default;

  ModeShapeBase(ModeShapeBase&&) noexcept = default;

  auto operator=(const ModeShapeBase&) -> ModeShapeBase& = default;

  auto operator=(ModeShapeBase&&) noexcept -> ModeShapeBase& = default;

  virtual auto buildModel(Array1D<Real> x_node_positions,
                          Array1D<Real> y_node_positions,
                          Array1D<Real> z_node_positions) -> void = 0;

  auto& gridModel() const { return _grid_model; }

  auto& gridModel() { return _grid_model; }

 protected:
  std::shared_ptr<GridModel> _grid_model;
};

template <unit::Length Unit>
class ModelShape : public ModeShapeBase {
 public:
  ModelShape(std::istream& model_info, std::istream& model_vertices,
             std::istream& model_elements) {
    _grid_model = std::make_shared<GridModel>();
    _grid_model->read(model_info, model_vertices, model_elements);
  }

  auto clone() const -> std::unique_ptr<Shape> override {
    return std::make_unique<ModelShape>(*this);
  }

  auto isInside(xfdtd::Real x, xfdtd::Real y, xfdtd::Real z,
                xfdtd::Real eps) const -> bool override {
    x = standardToUnit(x);
    y = standardToUnit(y);
    z = standardToUnit(z);
    eps = standardToUnit(eps);
    return _grid_model->isInside(x, y, z, eps);
  }

  auto isInside(const xfdtd::Vector& vector,
                xfdtd::Real eps) const -> bool override {
    return isInside(vector.x(), vector.y(), vector.z(), eps);
  }

  auto wrappedCube() const -> std::unique_ptr<xfdtd::Cube> override {
    auto&& triangular_model_info = _grid_model->triangularModelInfo();
    return std::make_unique<xfdtd::Cube>(
        xfdtd::Vector{unitToStandard(triangular_model_info.minX()),
                      unitToStandard(triangular_model_info.minY()),
                      unitToStandard(triangular_model_info.minZ())},
        xfdtd::Vector{unitToStandard(triangular_model_info.sizeX()),
                      unitToStandard(triangular_model_info.sizeY()),
                      unitToStandard(triangular_model_info.sizeZ())});
  }

  auto buildModel(Array1D<Real> x_node_positions,
                  Array1D<Real> y_node_positions,
                  Array1D<Real> z_node_positions) -> void override {
    for (auto&& x : x_node_positions) {
      x = standardToUnit(x);
    }
    for (auto&& y : y_node_positions) {
      y = standardToUnit(y);
    }
    for (auto&& z : z_node_positions) {
      z = standardToUnit(z);
    }
    _grid_model->buildModel(std::move(x_node_positions),
                            std::move(y_node_positions),
                            std::move(z_node_positions));
  }

  static constexpr auto unitToStandard(xfdtd::Real value) -> xfdtd::Real {
    if constexpr (Unit == unit::Length::Meter) {
      return value;
    } else if constexpr (Unit == unit::Length::Millimeter) {
      return value * 1e-3;
    } else if constexpr (Unit == unit::Length::Micrometer) {
      return value * 1e-6;
    } else if constexpr (Unit == unit::Length::Nanometer) {
      return value * 1e-9;
    } else {
      static_assert(Unit == unit::Length::Meter, "Invalid unit::Length");
    }
  }

  static constexpr auto standardToUnit(xfdtd::Real value) -> xfdtd::Real {
    if constexpr (Unit == unit::Length::Meter) {
      return value;
    } else if constexpr (Unit == unit::Length::Millimeter) {
      return value * 1e3;
    } else if constexpr (Unit == unit::Length::Micrometer) {
      return value * 1e6;
    } else if constexpr (Unit == unit::Length::Nanometer) {
      return value * 1e9;
    } else {
      static_assert(Unit == unit::Length::Meter, "Invalid unit::Length");
    }
  }
};

}  // namespace xfdtd::model

#endif  // __XFDTD_MODEL_MODEL_SHAPE_HPP__
