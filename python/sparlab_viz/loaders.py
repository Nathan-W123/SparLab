"""Readers for the artefacts a SparLab run writes.

Nothing here recomputes physics. Every loader is a thin, validated reader, so a
missing or malformed file produces a clear error instead of a silently empty
plot.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd


class ResultError(RuntimeError):
    """Raised when a result directory is missing or inconsistent."""


def load_json(path: str) -> Dict[str, Any]:
    """Read a JSON document written by SparLab."""
    if not os.path.isfile(path):
        raise ResultError(f"{path} does not exist; run the corresponding case first")
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_csv(path: str, required: bool = True) -> Optional[pd.DataFrame]:
    """Read one of the CSV tables. Returns None when optional and absent."""
    if not os.path.isfile(path):
        if required:
            raise ResultError(f"{path} does not exist; run the corresponding case first")
        return None
    frame = pd.read_csv(path)
    if frame.empty:
        raise ResultError(f"{path} contains no rows")
    return frame


@dataclass
class Mesh:
    """Nodes, connectivity, prescribed DOFs and applied loads of one case."""

    nodes: np.ndarray          # (num_nodes, 2) [m]
    elements: np.ndarray       # (num_elements, nodes_per_element)
    element_type: str
    prescribed: List[Dict[str, Any]] = field(default_factory=list)
    load_cases: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def num_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def num_elements(self) -> int:
        return int(self.elements.shape[0])

    @property
    def extent(self) -> tuple:
        """(xmin, xmax, ymin, ymax) of the nodal coordinates [m]."""
        return (
            float(self.nodes[:, 0].min()),
            float(self.nodes[:, 0].max()),
            float(self.nodes[:, 1].min()),
            float(self.nodes[:, 1].max()),
        )

    @property
    def polygons(self) -> np.ndarray:
        """(num_elements, nodes_per_element, 2) array of element corner points."""
        return self.nodes[self.elements]

    @property
    def element_centroids(self) -> np.ndarray:
        return self.polygons.mean(axis=1)

    def constrained_nodes(self, component: Optional[str] = None) -> np.ndarray:
        """Indices of nodes with a prescribed DOF, optionally for one component."""
        ids = [
            entry["node"]
            for entry in self.prescribed
            if component is None or entry["component"] == component
        ]
        return np.unique(np.asarray(ids, dtype=int)) if ids else np.empty(0, dtype=int)

    def nodal_forces(self, load_case: str) -> np.ndarray:
        """(k, 3) array of [node, fx, fy] for the named load case [N]."""
        for case in self.load_cases:
            if case["name"] != load_case:
                continue
            rows = [
                (entry["node"], entry["fx_N"], entry["fy_N"])
                for entry in case["nodal_forces"]
            ]
            return np.asarray(rows, dtype=float) if rows else np.empty((0, 3))
        raise ResultError(f"load case '{load_case}' is not present in mesh.json")

    @property
    def load_case_names(self) -> List[str]:
        return [case["name"] for case in self.load_cases]


def load_mesh(directory: str) -> Mesh:
    """Read `mesh.json` from a result directory."""
    doc = load_json(os.path.join(directory, "mesh.json"))
    nodes = np.asarray(doc["nodes_m"], dtype=float)
    elements = np.asarray(doc["elements"], dtype=int)
    if nodes.ndim != 2 or nodes.shape[1] != 2:
        raise ResultError("mesh.json: nodes_m must be an array of [x, y] pairs")
    if elements.ndim != 2:
        raise ResultError("mesh.json: elements must be a rectangular array")
    if elements.size and (elements.min() < 0 or elements.max() >= nodes.shape[0]):
        raise ResultError("mesh.json: connectivity references a node outside the mesh")
    return Mesh(
        nodes=nodes,
        elements=elements,
        element_type=doc.get("element_type", "Quad4"),
        prescribed=doc.get("prescribed_dofs", []),
        load_cases=doc.get("load_cases", []),
    )


def load_summary(directory: str) -> Dict[str, Any]:
    """Read `summary.json` from a result directory."""
    return load_json(os.path.join(directory, "summary.json"))


@dataclass
class CaseResults:
    """Everything one result directory holds, loaded lazily where it is large."""

    directory: str
    mesh: Mesh
    summary: Dict[str, Any]

    def path(self, name: str) -> str:
        return os.path.join(self.directory, name)

    def has(self, name: str) -> bool:
        return os.path.isfile(self.path(name))

    # -- static fields -----------------------------------------------------
    def displacement(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"displacement_{_safe(load_case)}.csv"))

    def stress(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"stress_{_safe(load_case)}.csv"))

    def reactions(self, load_case: str) -> pd.DataFrame:
        return load_csv(self.path(f"reactions_{_safe(load_case)}.csv"))

    # -- modal -------------------------------------------------------------
    def modes(self, tag: str = "") -> Optional[pd.DataFrame]:
        suffix = f"_{_safe(tag)}" if tag else ""
        return load_csv(self.path(f"modes{suffix}.csv"), required=False)

    def mode_shapes(self, tag: str = "") -> Optional[pd.DataFrame]:
        suffix = f"_{_safe(tag)}" if tag else ""
        return load_csv(self.path(f"mode_shapes{suffix}.csv"), required=False)

    # -- topology ----------------------------------------------------------
    def history(self) -> Optional[pd.DataFrame]:
        return load_csv(self.path("history.csv"), required=False)

    def density(self) -> Optional[pd.DataFrame]:
        return load_csv(self.path("density_final.csv"), required=False)

    def density_history(self) -> Optional[pd.DataFrame]:
        return load_csv(self.path("density_history.csv"), required=False)

    @property
    def load_case_names(self) -> List[str]:
        return self.mesh.load_case_names

    @property
    def name(self) -> str:
        return self.summary.get("provenance", {}).get("case", os.path.basename(self.directory))

    @property
    def is_topology_run(self) -> bool:
        return "optimization_result" in self.summary


def load_case(directory: str) -> CaseResults:
    """Load the mesh and summary of a result directory."""
    if not os.path.isdir(directory):
        raise ResultError(f"{directory} is not a directory")
    return CaseResults(
        directory=directory,
        mesh=load_mesh(directory),
        summary=load_summary(directory),
    )


def _safe(name: str) -> str:
    """Mirror the file-name sanitisation done by the C++ ResultWriter."""
    return "".join(
        c if (c.isalnum() or c in "-_") else "_" for c in name
    ) or "unnamed"
