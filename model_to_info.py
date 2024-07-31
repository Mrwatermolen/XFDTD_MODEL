import meshio
import numpy as np
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert mesh file to txt file.")
    parser.add_argument("--mesh_file", type=str, help="The mesh file.")
    parser.add_argument("--save_path", default="./", type=str, help="The path to save the txt files.")
    
    args = parser.parse_args()
    
    mesh_file = args.mesh_file
    save_path = args.save_path
    if os.path.exists(mesh_file) == False:
        raise ValueError(f"The mesh file {mesh_file} does not exist.")
    mesh = meshio.read(mesh_file)

    # 节点坐标提取
    nodes = mesh.points
    Nn = nodes.shape[0]

    # 三角面元提取
    ele_triangle = mesh.cells_dict["triangle"]
    Ne = ele_triangle.shape[0]

    if os.path.exists(args.save_path) == False:
        print(f"Create directory {args.save_path} to save the txt files.")
        os.makedirs(args.save_path)
        
    np.savetxt(os.path.join(args.save_path, 'nodes.txt'), nodes, fmt = '%20.15f')
    np.savetxt(os.path.join(args.save_path, 'elements.txt'), ele_triangle, fmt = '%5d')
    np.savetxt(os.path.join(args.save_path, 'info.txt'), [Nn, Ne], fmt = "%d")
    