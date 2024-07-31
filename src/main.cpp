#include <iomanip>
#include <ios>
#include <random>
#include <sstream>

#include "argparse/argparse.hpp"
#include "fdtd_model.hpp"

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

auto main(int argc, char** argv) -> int {
  argparse::ArgumentParser program("fdtd_model_msh");
  program.add_argument("--delta_l", "-l")
      .help("the resolution of grid (m)")
      .required()
      .scan<'f', float>();
  program.add_argument("--model_info_dir", "-d")
      .help("model_info_dir")
      .required();
  program.add_argument("--save_dir", "-s").help("save_dir").default_value(".");
  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << "\n";
    std::cerr << program;
    return 1;
  }

  auto delta_l = program.get<float>("--delta_l");
  auto model_info_dir = program.get<std::string>("--model_info_dir");
  auto save_dir = program.get<std::string>("--save_dir");

  // print info
  {
    std::stringstream ss;
    ss << "Run with:\n";
    ss << "The resolution of grid: " << delta_l << "\n";
    ss << "Model info dir: " << model_info_dir << "\n";
    ss << "Save dir: " << save_dir << "\n";
    std::cout << ss.str();
  }

  std::chrono::steady_clock::time_point begin =
      std::chrono::steady_clock::now();
  auto model = FDTDModel{delta_l, delta_l, delta_l};
  model.read(model_info_dir);
  std::cout << model.toString() << "\n";

  std::stringstream ss;
  ss << randomString(6) << "_" << std::setprecision(4) << std::fixed
     << (delta_l * 1e3) << "mm" << ".msh";
  auto save_name = ss.str();

  std::filesystem::path save_path = save_dir;
  if (!std::filesystem::exists(save_path)) {
    std::cout << "Create save directory: " << save_path.string() << "\n";
    std::filesystem::create_directories(save_path);
  }
  save_path /= save_name;
  std::cout << "Save to: " << save_path.string() << "\n";
  model.saveToMsh(save_path.string());
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Elapsed time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end -
                                                                     begin)
                   .count()
            << "[ms]" << "\n";
  return 0;
}
