#include <filesystem>
#include <fstream>
#include <memory>

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
#include "xfdtd_model/model_object.hpp"
#include "xfdtd_model/model_shape.hpp"

auto modelScatter() -> void {
  using namespace std::string_view_literals;
  constexpr auto model_info_dir = "./tmp/model/info/HTV2/metal"sv;
  constexpr auto data_path_str = "./tmp/data/HTV2"sv;
  const auto data_path = std::filesystem::path{data_path_str};
  constexpr xfdtd::Real delta_l = 25e-3;

  using ModelShape = xfdtd::model::ModelShape<xfdtd::unit::Length::Millimeter>;
  // try to open info.txt, nodes.txt, and elements.txt
  auto model_info_dir_path = std::filesystem::path{model_info_dir};
  auto info_path = model_info_dir_path / "info.txt";
  auto nodes_path = model_info_dir_path / "nodes.txt";
  auto elements_path = model_info_dir_path / "elements.txt";
  auto info_file = std::ifstream{info_path};
  auto nodes_file = std::ifstream{nodes_path};
  auto elements_file = std::ifstream{elements_path};
  if (!info_file || !nodes_file || !elements_file) {
    throw std::runtime_error("Failed to open model info files");
  }

  auto model_shape =
      std::make_unique<ModelShape>(info_file, nodes_file, elements_file);
  auto wrapped_cube = model_shape->wrappedCube();
  std::cout << "Model Shape: " << wrapped_cube->toString() << "\n";
  auto domain_shape = std::make_unique<xfdtd::Cube>(
      wrapped_cube->origin() -
          xfdtd::Vector{10 * delta_l, 10 * delta_l, 10 * delta_l},
      wrapped_cube->size() +
          xfdtd::Vector{20 * delta_l, 20 * delta_l, 20 * delta_l});

  auto domain = std::make_shared<xfdtd::Object>(
      "domain", domain_shape->wrappedCube(), xfdtd::Material::createAir());

  auto model = std::make_shared<xfdtd::model::ModelObject>(
      "model", std::move(model_shape), xfdtd::Material::createPec());

  constexpr auto l_min{delta_l * 20};
  constexpr auto f_max{3e8 / l_min};  // max frequency: 5 GHz in dl = 3e-3
  constexpr auto tau{l_min / 6e8};
  constexpr auto t_0{4.5 * tau};
  constexpr xfdtd::Index tfsf_start{11};
  auto tfsf{std::make_shared<xfdtd::TFSF3D>(
      tfsf_start, tfsf_start, tfsf_start, xfdtd::constant::PI / 2, 0, 1,
      xfdtd::Waveform::gaussian(tau, t_0))};

  auto movie_ez_xy{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(
              xfdtd::Vector{domain_shape->originX(),
                            wrapped_cube->originY() + 3 * delta_l,
                            domain_shape->originZ()},
              xfdtd::Vector{domain_shape->sizeX(), delta_l,
                            domain_shape->sizeZ()}),
          xfdtd::EMF::Field::EZ, "", ""),
      15, "movie_ez_xz", (data_path / "movie_ez_xz").string())};

  auto movie_ez_mid_xy{std::make_shared<xfdtd::MovieMonitor>(
      std::make_unique<xfdtd::FieldMonitor>(
          std::make_unique<xfdtd::Cube>(
              xfdtd::Vector{domain_shape->originX(),
                            wrapped_cube->originY() +
                                wrapped_cube->sizeY() / 2 + 5 * delta_l,
                            domain_shape->originZ()},
              xfdtd::Vector{domain_shape->sizeX(), delta_l,
                            domain_shape->sizeZ()}),
          xfdtd::EMF::Field::EZ, "", ""),
      15, "movie_ez_mid_xz", (data_path / "movie_ez_mid_xz").string())};

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
  s.addMonitor(movie_ez_xy);
  s.addMonitor(movie_ez_mid_xy);
  s.run(800);

  // Save to msh
  if (!xfdtd::MpiSupport::instance().isRoot()) {
    return;
  }

  // {
  //   std::stringstream ss;
  //   ss << "Input y to save to msh: ";
  //   std::string input;
  //   std::cout << ss.str();

  //   std::cin >> input;
  //   if (input != "y") {
  //     return;
  //   }

  //   auto name = model_info_dir.substr(model_info_dir.find_last_of('/') + 1);
  //   ss.str("");
  //   ss << data_path_str << "/" << name << "_" << std::setprecision(4)
  //      << std::fixed << (delta_l * 1e3) << "mm" << ".msh";
  //   auto save_path = ss.str();
  //   std::cout << "Save to: " << save_path << "\n";
  //   auto model_shape =
  //       dynamic_cast<ModelShape<LengthUnit::Millimeter>*>(model->shape().get());
  //   if (!model_shape) {
  //     throw std::runtime_error("Failed to cast model shape");
  //   }

  //   model_shape->model()->saveToMsh(save_path);
  // }
}

int main() { modelScatter(); }
