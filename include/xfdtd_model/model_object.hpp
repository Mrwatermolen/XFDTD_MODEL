#ifndef __XFDTD_MODEL_MODEL_OBJECT_HPP__
#define __XFDTD_MODEL_MODEL_OBJECT_HPP__

#include <xfdtd/object/object.h>

#include <string_view>

#include "xfdtd_model/exception.h"
#include "xfdtd_model/model_shape.hpp"

namespace xfdtd::model {

class XFDTDModelModelObjectException : public XFDTDModelException {
 public:
  using XFDTDModelException::XFDTDModelException;
};

class ModelObject : public Object {
 public:
  ModelObject(std::string_view name, std::unique_ptr<ModeShapeBase> shape,
              const std::shared_ptr<Material>& material, Real factor = 1.0)
      : Object{name.data(), std::move(shape), material},
        _factor{factor},
        _model_shape{dynamic_cast<ModeShapeBase*>(shapePtr())} {
    if (_model_shape == nullptr) {
      throw XFDTDModelModelObjectException{
          "Failed to cast shape to ModeShapeBase"};
    }
  }

  ModelObject(const ModelObject&) = delete;

  ModelObject(ModelObject&&) noexcept = default;

  ModelObject& operator=(const ModelObject&) = delete;

  ModelObject& operator=(ModelObject&&) noexcept = default;

  ~ModelObject() override = default;

  // I hate <const GridSpace>!!!
  auto init(std::shared_ptr<const GridSpace> grid_space,
            std::shared_ptr<CalculationParam> calculation_param,
            std::shared_ptr<EMF> emf) -> void override {
    Object::init(grid_space, calculation_param, emf);
    _model_shape->buildModel(grid_space->eNodeX(), grid_space->eNodeY(),
                             grid_space->eNodeZ(), _factor);
  }

  auto& modelShape() const { return *_model_shape; }

  /**
   * @brief The factor for scaling the resolution of the model
   * 
   * @return Real
   */
  auto factor() const { return _factor; }

  auto setFactor(Real factor) { _factor = factor; }

 private:
  ModeShapeBase* _model_shape{};
  Real _factor{1.0}; // For scaling the resolution of the model
};

}  // namespace xfdtd::model

#endif  // __XFDTD_MODEL_MODEL_OBJECT_HPP__
