import sys
import numpy as np

def main():
    if len(sys.argv) < 2:
        raise ValueError("Usage: python script.py <mesh_dir>")

    mesh_dir = sys.argv[1]
    if not mesh_dir.endswith("/"):
        mesh_dir += "/"

    nodhn_file = mesh_dir + "nodhn.out"
    nodhn = np.loadtxt(nodhn_file)

    min_depth = -20.0
    depth = (-1) * nodhn
    depth = np.where(depth > min_depth, min_depth, depth)

    with open("depth@node.out", "w") as f:
        np.savetxt(f, depth, fmt="%f")

if __name__ == "__main__":
    main()

