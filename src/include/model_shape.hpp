
#include <memory>
#include <string>
#include <string_view>

#include "fdtd_model.hpp"
#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/shape/shape.h"

enum class LengthUnit { Meter, Millimeter, Micrometer, Nanometer };

template <LengthUnit Unit>
class ModelShape : public xfdtd::Shape {
 public:
  ModelShape(std::string_view model_info_path, xfdtd::Real dl)
      : _model_info_path{model_info_path},
        _model{std::make_shared<FDTDModel>(dl, dl, dl)} {
    _model->read(model_info_path);
  }

  ~ModelShape() override = default;

  auto clone() const -> std::unique_ptr<xfdtd::Shape> override {
    return std::make_unique<ModelShape>(*this);
  }

  auto isInside(xfdtd::Real x, xfdtd::Real y, xfdtd::Real z,
                xfdtd::Real eps) const -> bool override {
    return _model->isInside(standardToUnit(x), standardToUnit(y),
                            standardToUnit(z), standardToUnit(eps));
  }

  auto isInside(const xfdtd::Vector& vector,
                xfdtd::Real eps) const -> bool override {
    return _model->isInside(standardToUnit(vector.x()),
                            standardToUnit(vector.y()),
                            standardToUnit(vector.z()), standardToUnit(eps));
  }

  auto wrappedCube() const -> std::unique_ptr<xfdtd::Cube> override {
    auto domain = _model->domain();

    return std::make_unique<xfdtd::Cube>(
        xfdtd::Vector{unitToStandard(domain.start().x()),
                      unitToStandard(domain.start().y()),
                      unitToStandard(domain.start().z())},
        xfdtd::Vector{unitToStandard(domain.size().x()),
                      unitToStandard(domain.size().y()),
                      unitToStandard(domain.size().z())});
  }

  auto releaseModel() -> void { _model.reset(); }

  auto model() const { return _model; }

  auto model() { return _model; }

  static constexpr auto unitToStandard(xfdtd::Real value) -> xfdtd::Real {
    if constexpr (Unit == LengthUnit::Meter) {
      return value;
    } else if constexpr (Unit == LengthUnit::Millimeter) {
      return value * 1e-3;
    } else if constexpr (Unit == LengthUnit::Micrometer) {
      return value * 1e-6;
    } else if constexpr (Unit == LengthUnit::Nanometer) {
      return value * 1e-9;
    } else {
      static_assert(Unit == LengthUnit::Meter, "Invalid LengthUnit");
    }
  }

  static constexpr auto standardToUnit(xfdtd::Real value) -> xfdtd::Real {
    if constexpr (Unit == LengthUnit::Meter) {
      return value;
    } else if constexpr (Unit == LengthUnit::Millimeter) {
      return value * 1e3;
    } else if constexpr (Unit == LengthUnit::Micrometer) {
      return value * 1e6;
    } else if constexpr (Unit == LengthUnit::Nanometer) {
      return value * 1e9;
    } else {
      static_assert(Unit == LengthUnit::Meter, "Invalid LengthUnit");
    }
  }

 private:
  std::string _model_info_path;
  std::shared_ptr<FDTDModel> _model;
};
