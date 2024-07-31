#include <filesystem>
#include <memory>

#include "model_shape.hpp"
#include "xfdtd/boundary/pml.h"
#include "xfdtd/common/constant.h"
#include "xfdtd/common/type_define.h"
#include "xfdtd/coordinate_system/coordinate_system.h"
#include "xfdtd/material/material.h"
#include "xfdtd/monitor/field_monitor.h"
#include "xfdtd/monitor/movie_monitor.h"
#include "xfdtd/object/object.h"
#include "xfdtd/parallel/mpi_support.h"
#include "xfdtd/shape/cube.h"
#include "xfdtd/simulation/simulation.h"
#include "xfdtd/waveform_source/tfsf_3d.h"

auto modelScatter() -> void {
  using namespace std::string_view_literals;
  constexpr auto model_info_dir =
      "./tmp/model/info/HTV2"sv;
  constexpr auto data_path_str =
      "./tmp/data/HTV2"sv;
  const auto data_path = std::filesystem::path{data_path_str};
  constexpr xfdtd::Real delta_l = 30e-3;

  auto model_shape = std::make_unique<ModelShape<LengthUnit::Millimeter>>(
      model_info_dir,
      ModelShape<LengthUnit::Millimeter>::standardToUnit(delta_l));
  auto wrapped_cube = model_shape->wrappedCube();
  std::cout << "Model Shape: " << wrapped_cube->toString() << "\n";
  auto domain_shape = std::make_unique<xfdtd::Cube>(
      wrapped_cube->origin() -
          xfdtd::Vector{10 * delta_l, 10 * delta_l, 10 * delta_l},
      wrapped_cube->size() +
          xfdtd::Vector{20 * delta_l, 20 * delta_l, 20 * delta_l});

  auto domain = std::make_shared<xfdtd::Object>(
      "domain", domain_shape->wrappedCube(), xfdtd::Material::createAir());

  auto model = std::make_shared<xfdtd::Object>("model", std::move(model_shape),
                                               xfdtd::Material::createPec());

  constexpr auto l_min{delta_l * 20};
  constexpr auto f_max{3e8 / l_min};  // max frequency: 5 GHz in dl = 3e-3
  constexpr auto tau{l_min / 6e8};
  constexpr auto t_0{4.5 * tau};
  constexpr std::size_t tfsf_start{static_cast<size_t>(11)};
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, xfdtd::constant::PI / 2, 0, 1,
      xfdtd::Waveform::gaussian(tau, t_0))};

  //   auto movie_ex_xz{std::make_shared<xfdtd::MovieMonitor>(
  //       std::make_unique<xfdtd::FieldMonitor>(
  //           std::make_unique<xfdtd::Cube>(
  //               xfdtd::Vector{domain_shape->originX(), 0,
  //                             domain_shape->originZ()},
  //               xfdtd::Vector{domain_shape->sizeX(), delta_l,
  //                             domain_shape->sizeZ()}),
  //           xfdtd::EMF::Field::EX, "", ""),
  //       15, "movie_ex_xz", (data_path / "movie_ex_xz").string())};

  auto movie_ez_xy{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(
              xfdtd::Vector{domain_shape->originX(), domain_shape->originY(),
                            0},
              xfdtd::Vector{domain_shape->sizeX(), domain_shape->sizeY(),
                            delta_l}),
          xfdtd::EMF::Field::EZ, "", ""),
      15, "movie_ez_xy", (data_path / "movie_ez_xy").string())};

  auto movie_ez_mid_xy{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(
              xfdtd::Vector{
                  domain_shape->originX(), domain_shape->originY(),
                  domain_shape->originZ() + domain_shape->sizeZ() / 2},
              xfdtd::Vector{domain_shape->sizeX(), domain_shape->sizeY(),
                            delta_l}),
          xfdtd::EMF::Field::EZ, "", ""),
      15, "movie_ez_mid_xy", (data_path / "movie_ez_mid_xy").string())};

  auto s{xfdtd::Simulation{delta_l, delta_l, delta_l, 0.9,
                           xfdtd::ThreadConfig{1, 1, 1}}};
  s.addObject(domain);
  s.addObject(model);
  s.addWaveformSource(tfsf);
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::XP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::YP));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZN));
  s.addBoundary(std::make_shared<xfdtd::PML>(8, xfdtd::Axis::Direction::ZP));
  //   s.addMonitor(movie_ex_xz);
  s.addMonitor(movie_ez_xy);
  s.addMonitor(movie_ez_mid_xy);
  s.run(2400);

  // Save to msh
  if (!xfdtd::MpiSupport::instance().isRoot()) {
    return;
  }

  {
    std::stringstream ss;
    ss << "Input y to save to msh: ";
    std::string input;
    std::cout << ss.str();

    std::cin >> input;
    if (input != "y") {
      return;
    }

    auto name = model_info_dir.substr(model_info_dir.find_last_of('/') + 1);
    ss.str("");
    ss << data_path_str << "/" << name << "_" << std::setprecision(4)
       << std::fixed << (delta_l * 1e3) << "mm" << ".msh";
    auto save_path = ss.str();
    std::cout << "Save to: " << save_path << "\n";
    auto model_shape =
        dynamic_cast<ModelShape<LengthUnit::Millimeter>*>(model->shape().get());
    if (!model_shape) {
      throw std::runtime_error("Failed to cast model shape");
    }

    model_shape->model()->saveToMsh(save_path);
  }
}

int main() { modelScatter(); }
