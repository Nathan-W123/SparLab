import json
import numpy as np


def write_bin(path, X, F):
    with open(path, "wb") as f:
        np.array([len(X), len(F)], dtype=np.uint32).tofile(f)
        np.asarray(X, dtype=np.float32).tofile(f)
        np.asarray(F, dtype=np.uint32).tofile(f)


def write_index(directory, count, scale, **extra):
    json.dump({"count": count, "scale": scale, **extra}, open(f"{directory}/index.json", "w"), indent=1)
