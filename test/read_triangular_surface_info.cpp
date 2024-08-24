#include <ase_reader/ase_reader.hpp>
#include <iostream>
#include <sstream>
#include <xfdtd_model/triangular_model_info.hpp>

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 1;
  }

  xfdtd::model::TriangularModelInfo info;
  ase_reader::ASEReader ase;
  ase.read(argv[1]);

  {
    std::stringstream info_ss;
    std::stringstream vertices_ss;
    std::stringstream faces_ss;

    ase.write(info_ss, vertices_ss, faces_ss);

    info.read(info_ss, vertices_ss, faces_ss);
  }

  std::cout << "Region: " << "minX=" << info.minX() << ", minY=" << info.minY()
            << ", minZ=" << info.minZ() << ", sizeX=" << info.sizeX()
            << ", sizeY=" << info.sizeY() << ", sizeZ=" << info.sizeZ()
            << ", maxX=" << info.maxX() << ", maxY=" << info.maxY()
            << ", maxZ=" << info.maxZ() << std::endl;

  std::cout << "Vertices: " << info.numVertices() << std::endl;
  std::cout << "Faces: " << info.numElements() << std::endl;

  return 0;
}
