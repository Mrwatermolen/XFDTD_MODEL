#include <argparse/argparse.hpp>
#include <ase_reader/ase_reader.hpp>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <xfdtd_model/grid_model.hpp>

static auto randomString(const auto& length) {
  static thread_local std::mt19937 generator(std::random_device{}());

  static auto rand_char = [&]() -> char {
    const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const std::size_t max_index = (sizeof(charset) - 1);
    std::uniform_int_distribution<std::size_t> distribution(0, max_index);
    return charset[distribution(generator) % max_index];
  };

  std::string str(length, 0);
  std::generate_n(str.begin(), length, rand_char);
  return str;
}

int main(int argc, char** argv) {
  argparse::ArgumentParser program("fdtd_model_msh");

  program.add_argument("--delta_l", "-l")
      .help("the resolution of grid")
      .required()
      .scan<'f', float>();

  program.add_argument("--model_info_dir", "-d").help("model_info_dir");

  program.add_argument("--ase_file", "-a").help("ase_file");

  program.add_argument("--save_dir", "-s").help("save_dir").default_value(".");

  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << "\n";
    std::cerr << program;
    return 1;
  }

  auto delta_l = program.get<float>("--delta_l");
  std::string ase_file{};
  std::string model_info_dir{};
  try {
    ase_file = program.get<std::string>("--ase_file");
  } catch (const std::logic_error& err) {
    try {
      model_info_dir = program.get<std::string>("--model_info_dir");
    } catch (const std::logic_error& err) {
      std::cerr << "Please provide either --ase_file or --model_info_dir\n";
      std::cerr << program;
      return 1;
    }
  }

  auto save_dir = program.get<std::string>("--save_dir");

  xfdtd::model::GridModel grid_model;

  {
    std::stringstream info_ss;
    std::stringstream vertices_ss;
    std::stringstream faces_ss;

    if (!ase_file.empty()) {
      ase_reader::ASEReader ase;
      ase.read(ase_file);
      ase.write(info_ss, vertices_ss, faces_ss);
    }

    if (!model_info_dir.empty()) {
      auto info_path = std::filesystem::path(model_info_dir) / "info.txt";
      auto nodes_path = std::filesystem::path(model_info_dir) / "nodes.txt";
      auto elements_path =
          std::filesystem::path(model_info_dir) / "elements.txt";
      if (!std::filesystem::exists(info_path)) {
        std::cerr << "info.txt not found\n";
        return 1;
      }

      if (!std::filesystem::exists(nodes_path)) {
        std::cerr << "nodes.txt not found\n";
        return 1;
      }

      if (!std::filesystem::exists(elements_path)) {
        std::cerr << "elements.txt not found\n";
        return 1;
      }

      auto is = std::ifstream(info_path);
      info_ss << is.rdbuf();
      is.close();
      is = std::ifstream(nodes_path);
      vertices_ss << is.rdbuf();
      is.close();
      is = std::ifstream(elements_path);
      faces_ss << is.rdbuf();
    }

    grid_model.read(info_ss, vertices_ss, faces_ss);
  }

  {
    auto info = grid_model.triangularModelInfo();
    std::cout << "Region: " << "minX=" << info.minX()
              << ", minY=" << info.minY() << ", minZ=" << info.minZ()
              << ", sizeX=" << info.sizeX() << ", sizeY=" << info.sizeY()
              << ", sizeZ=" << info.sizeZ() << ", maxX=" << info.maxX()
              << ", maxY=" << info.maxY() << ", maxZ=" << info.maxZ()
              << std::endl;

    std::cout << "Vertices: " << info.numVertices() << std::endl;
    std::cout << "Faces: " << info.numElements() << std::endl;
  }

  {
    std::cout << "Build grid model with delta_l=" << delta_l << std::endl;
    grid_model.buildUniform(delta_l, delta_l, delta_l);
    std::stringstream ss;
    ss << randomString(6) << "_" << std::setprecision(4) << std::fixed
       << (delta_l) << ".msh";
    auto save_name = ss.str();
    auto save_path = std::filesystem::path(save_dir) / save_name;
    if (!std::filesystem::exists(save_path.parent_path())) {
      std::cout << "Create save directory: " << save_path.parent_path().string()
                << "\n";
      std::filesystem::create_directories(save_path.parent_path());
    }

    std::cout << "Save to: " << save_path.string() << "\n";
    std::ofstream os(save_path);

    grid_model.writeToMsh(os);
  }
}
