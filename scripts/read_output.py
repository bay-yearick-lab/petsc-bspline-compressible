"""Load a field written by the channel-flow driver (<var>_<step>.out).

Each file holds one line per B-spline coefficient: x, y, z, value, ordered with
x fastest, then y, then z. Returns arrays shaped (nz, ny, nx).

Usage: python scripts/read_output.py u_0000000100.out [--plot slice.png]
"""
import argparse

import numpy as np


def load_field(path):
    """Return (x, y, z, f) with f shaped (nz, ny, nx)."""
    data = np.loadtxt(path)
    x = np.unique(data[:, 0])
    y = np.unique(data[:, 1])
    z = np.unique(data[:, 2])
    f = data[:, 3].reshape(len(z), len(y), len(x))
    return x, y, z, f


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("path")
    ap.add_argument("--plot", help="write a streamwise/wall-normal slice to this PNG")
    args = ap.parse_args()

    x, y, z, f = load_field(args.path)
    print(f"{args.path}: nx={len(x)} ny={len(y)} nz={len(z)} "
          f"min={f.min():.6e} max={f.max():.6e} mean={f.mean():.6e}")

    if args.plot:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(10, 3))
        im = ax.pcolormesh(x, y, f[len(z) // 2], shading="auto", cmap="viridis")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im, ax=ax)
        fig.savefig(args.plot, dpi=200, bbox_inches="tight")


if __name__ == "__main__":
    main()
