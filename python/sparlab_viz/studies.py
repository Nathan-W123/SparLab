"""Figures for the verification studies, the design study and the benchmark.

These read the CSV / JSON files written by `sparlab_verify`, `sparlab_bench` and
the study script. Like `plots`, nothing here recomputes physics; the only derived
values are ratios and least-squares slopes, and each is labelled as such.
"""

from __future__ import annotations

import glob
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from . import style as st
from .loaders import load_csv, load_json


# ---------------------------------------------------------------------------
# Verification: sensitivity finite differences
# ---------------------------------------------------------------------------
def plot_sensitivity_check(directory: str, path: str) -> str:
    """Analytical vs finite-difference gradients, and error vs step size."""
    steps = load_csv(os.path.join(directory, "sensitivity_steps.csv"))
    elements = load_csv(os.path.join(directory, "sensitivity_elements.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get("sensitivity", {})
    tolerance = summary.get("tolerance", 1e-5)
    best_step = summary.get("best_step")

    fig, axes = st.figure(7.4, 6.4, nrows=2, ncols=1)

    # --- panel 1: error vs step size (3 series: within the all-pairs cap) ---
    ax = axes[0]
    st.require_scatter_series(3, "error measures")
    included = steps["num_excluded"].iloc[0]
    for slot, (column, label) in enumerate(
        [
            ("max_relative_error[-]", "max relative error over elements"),
            ("rms_relative_error[-]", "RMS relative error over elements"),
            ("directional_relative_error[-]", "relative error of the directional derivative"),
        ]
    ):
        ax.plot(steps["step[-]"], steps[column], "o-", color=st.series_color(slot),
                label=label)
    ax.axhline(tolerance, color=st.INK_MUTED, linewidth=1.0, linestyle="--")
    ax.text(
        steps["step[-]"].min(), tolerance, f" pass tolerance {tolerance:g}",
        fontsize=7.8, color=st.INK_SECONDARY, va="bottom",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("central-difference step h on the design variable [-]")
    ax.set_ylabel("relative error [-]")
    st.title(
        ax, "Topology sensitivity: analytical gradient vs central differences",
        f"{int(steps['num_tested'].iloc[0])} elements tested, {int(included)} "
        "excluded at active bounds; the error floor is where round-off in the "
        "difference of two compliances overtakes the truncation error",
    )
    st.legend(ax, loc="upper left")

    # --- panel 2: element-by-element comparison at the best step ------------
    ax = axes[1]
    subset = elements
    if best_step is not None:
        subset = elements[np.isclose(elements["step[-]"], best_step)]
    tested = subset[subset["excluded[-]"] < 0.5]
    excluded = subset[subset["excluded[-]"] > 0.5]

    ax.plot(tested["analytical[J]"], tested["finite_difference[J]"], "o",
            color=st.series_color(0), markersize=4.5,
            label=f"tested elements ({len(tested)})")
    if not excluded.empty:
        ax.plot(excluded["analytical[J]"],
                np.zeros(len(excluded)), "x", color=st.series_color(1),
                markersize=5, label=f"excluded at bounds ({len(excluded)})")
    limits = [
        float(min(tested["analytical[J]"].min(), tested["finite_difference[J]"].min())),
        float(max(tested["analytical[J]"].max(), tested["finite_difference[J]"].max())),
    ]
    span = limits[1] - limits[0]
    limits = [limits[0] - 0.05 * span, limits[1] + 0.05 * span]
    ax.plot(limits, limits, "-", color=st.INK_MUTED, linewidth=1.0,
            label="exact agreement")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    st.limit_ticks(ax, x=4, y=5)
    ax.set_xlabel(r"analytical $\partial c/\partial x_e$ [J]")
    ax.set_ylabel(r"central difference $\partial c/\partial x_e$ [J]")
    worst = float(tested["relative_error[-]"].max()) if not tested.empty else float("nan")
    st.title(
        ax, f"element-by-element agreement at h = {best_step:g}"
        if best_step else "element-by-element agreement",
        f"worst relative error {worst:.3e}",
    )
    st.legend(ax, loc="lower right")

    st.annotate_note(
        fig,
        "Excluded elements are passive variables (lower bound equals upper "
        "bound) or free variables within h of 0 or 1, where a central "
        "difference would leave the feasible box. They are plotted at zero on "
        "the vertical axis only to show how many there are.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: mesh convergence and beam theory
# ---------------------------------------------------------------------------
def plot_mesh_convergence(directory: str, path: str) -> str:
    """Tip deflection and discretisation error vs element size."""
    table = load_csv(os.path.join(directory, "mesh_convergence.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("mesh_convergence", {})
    orders = block.get("observed_convergence_order_tip_deflection", {})

    poissons = sorted(table["poisson[-]"].unique())
    fig, axes = st.figure(7.4, 6.6, nrows=2, ncols=1)

    # --- panel 1: tip deflection vs element size, with both references -----
    ax = axes[0]
    for slot, nu in enumerate(poissons):
        sub = table[table["poisson[-]"] == nu].sort_values("h[m]")
        ax.plot(sub["h[m]"], np.abs(sub["tip_mean[m]"]), "o-",
                color=st.series_color(slot), label=f"FEM, nu = {nu:g}")
        ax.axhline(abs(float(sub["timoshenko[m]"].iloc[0])),
                   color=st.series_color(slot), linewidth=0.9, linestyle=":")
    ax.axhline(abs(float(table["euler_bernoulli[m]"].iloc[0])),
               color=st.INK_MUTED, linewidth=1.0, linestyle="--",
               label="Euler-Bernoulli (bending only)")
    ax.text(
        float(table["h[m]"].min()), abs(float(table["timoshenko[m]"].iloc[0])),
        " Timoshenko (bending + shear), dotted per nu", fontsize=7.6,
        color=st.INK_SECONDARY, va="bottom",
    )
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("|tip deflection| [m]")
    st.title(
        ax, "Cantilever tip deflection vs mesh size",
        "L = 1 m, h = 0.1 m (L/h = 10), t = 10 mm, E = 70 GPa, 1 kN tip "
        "resultant; deflection is the mean u_y over the tip edge",
    )
    st.legend(ax, loc="lower left")

    # --- panel 2: self-convergence error, i.e. verification ----------------
    ax = axes[1]
    for slot, nu in enumerate(poissons):
        sub = table[table["poisson[-]"] == nu].sort_values("h[m]", ascending=False)
        finest = float(sub["tip_mean[m]"].iloc[-1])
        error = np.abs(sub["tip_mean[m]"] - finest) / abs(finest)
        mask = error > 0
        ax.plot(sub["h[m]"][mask], error[mask], "o-", color=st.series_color(slot),
                label=f"nu = {nu:g} (observed order "
                      f"{orders.get(f'{nu:g}', float('nan')):.2f})")
    # Second-order reference slope.
    h_ref = np.array([float(table["h[m]"].min()) * 2.0, float(table["h[m]"].max())])
    anchor = 1e-3
    ax.plot(h_ref, anchor * (h_ref / h_ref[0]) ** 2, "-", color=st.INK_MUTED,
            linewidth=1.0, label=r"reference slope $h^2$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("relative error vs the finest mesh [-]")
    st.title(
        ax, "Self-convergence of the discretisation error",
        "measured against the finest mesh, which isolates discretisation error "
        "from the beam-theory modelling gap",
    )
    st.legend(ax, loc="lower left")

    st.annotate_note(
        fig,
        "Comparing against the finest mesh is verification (does the "
        "discretisation converge, and at what order); comparing against beam "
        "theory in the upper panel is validation (are the modelling assumptions "
        "appropriate), where a finite gap is expected because 1-D beam theory "
        "omits the Poisson coupling that plane stress retains.",
    )
    return st.save_figure(fig, path)


def plot_modal_convergence(directory: str, path: str) -> str:
    """Computed cantilever frequencies against Euler-Bernoulli theory."""
    table = load_csv(os.path.join(directory, "modal_convergence.csv"))
    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.5)

    ax = axes[0]
    series = [
        ("f1_fem[Hz]", "f1_theory[Hz]", "o", "bending mode 1"),
        ("f2_fem[Hz]", "f2_theory[Hz]", "s", "bending mode 2"),
        ("f_axial_fem[Hz]", "f_axial_theory[Hz]", "^", "axial mode 1"),
    ]
    for slot, (column, reference, marker, label) in enumerate(series):
        if column not in table:
            continue
        ax.plot(table["num_dofs"], table[column], marker + "-",
                color=st.series_color(slot), label=f"{label}, FEM")
        ax.axhline(float(table[reference].iloc[0]), color=st.series_color(slot),
                   linewidth=0.9, linestyle="--", label=f"{label}, theory")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("frequency [Hz]")
    st.title(
        ax, "Cantilever natural frequencies vs mesh size",
        "L = 1 m, h = 0.1 m, t = 10 mm, E = 70 GPa, rho = 2700 kg/m^3, nu = 0; "
        "consistent mass matrix. Bending references are Euler-Bernoulli; the "
        "axial reference is fixed-free rod theory",
    )
    st.legend(ax, loc="center right", ncol=2)

    ax = axes[1]
    errors = [
        ("f1_rel_error[-]", "o", "bending mode 1 vs Euler-Bernoulli"),
        ("f2_rel_error[-]", "s", "bending mode 2 vs Euler-Bernoulli"),
        ("f_axial_rel_error[-]", "^", "axial mode 1 vs rod theory"),
    ]
    for slot, (column, marker, label) in enumerate(errors):
        if column not in table:
            continue
        ax.plot(table["num_dofs"], table[column], marker + "-",
                color=st.series_color(slot), label=label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("relative difference [-]")
    st.title(
        ax, "difference from the analytical references",
        "the axial mode agrees to round-off because plane stress at nu = 0 IS a "
        "uniaxial bar; the bending modes sit a little below Euler-Bernoulli "
        "because the 2-D model includes shear deformation and rotary inertia "
        "that beam theory omits, and the gap grows with mode number",
    )
    st.legend(ax, loc="center right")
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Runtime scaling
# ---------------------------------------------------------------------------
def plot_runtime_scaling(directory: str, path: str, stem: str = "runtime_scaling",
                         dim: int = 2) -> str:
    """Wall-clock cost of each phase against problem size.

    `stem` selects the benchmark files (`runtime_scaling` for the Q4 plate,
    `runtime_scaling_3d` for the Hex8 block) and `dim` the wording.
    """
    table = load_csv(os.path.join(directory, f"{stem}.csv"))
    summary = load_json(os.path.join(directory, f"{stem}.json"))
    exponents = summary.get("scaling_exponent_vs_dofs", {})
    element = "Hex8" if dim == 3 else "Q4"

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.6)

    ax = axes[0]
    phases = [
        ("assemble[s]", "assemble K", "assemble"),
        ("factorize[s]", "sparse Cholesky factorisation", "factorize"),
        ("solve[s]", "one back-substitution", "solve"),
        ("objective_gradient[s]", "objective + gradient (one optimiser iteration)",
         "objective_gradient"),
    ]
    for slot, (column, label, key) in enumerate(phases):
        exponent = exponents.get(key)
        suffix = f"  (slope {exponent:.2f})" if exponent is not None else ""
        ax.plot(table["num_dofs"], table[column], "o-", color=st.series_color(slot),
                label=label + suffix)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel("wall-clock time [s]")
    st.title(
        ax, f"Runtime scaling of the solver phases ({element})",
        f"{summary.get('compiler', 'unknown compiler')}, "
        f"{summary.get('build_type', 'unknown build')}, "
        f"{summary.get('linear_solver', 'direct solver')}; "
        f"{summary.get('estimator', 'minimum over repeats')}",
    )
    st.legend(ax, loc="upper left")

    ax = axes[1]
    ax.plot(table["num_dofs"], table["stiffness_nonzeros"] / table["num_dofs"],
            "o-", color=st.series_color(2), label="stored entries per free DOF")
    ax.set_xscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("nonzeros / DOF [-]")
    st.title(
        ax, "sparsity of the reduced stiffness matrix",
        f"flat with size, as expected for a fixed-stencil structured {element} mesh",
    )
    st.legend(ax, loc="lower right")

    if dim == 3:
        asymptote = ("a 3-D sparse Cholesky with a fill-reducing ordering is close "
                     "to O(n^2) asymptotically (nested dissection), which is what "
                     "limits the direct solver on solid problems; the multigrid "
                     "solver compared in solver_scaling.png removes that limit")
    else:
        asymptote = ("a 2-D sparse Cholesky with a fill-reducing ordering is close "
                     "to O(n^1.5) asymptotically")
    st.annotate_note(
        fig,
        "Slopes are least-squares fits of log(time) against log(DOFs) over the "
        f"largest three sizes. Assembly is O(n); {asymptote}, and at these sizes "
        "the measured slope also carries cache effects.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Linear solvers: sparse Cholesky vs multigrid CG vs Jacobi CG
# ---------------------------------------------------------------------------
#: Colour slot per linear solver. The solver is the entity the colour names, so
#: it keeps its slot in every panel of every solver figure.
SOLVER_SLOTS = {"ldlt": 0, "amg": 1, "jacobi": 2}
SOLVER_NAMES = {
    "ldlt": "sparse Cholesky (LDL^T)",
    "amg": "multigrid-preconditioned CG",
    "jacobi": "Jacobi-preconditioned CG",
}
ELEMENT_MARKERS = {"Quad4": "o", "Tri3": "v", "Hex8": "s", "Tet4": "^", "Tet10": "D"}
ELEMENT_NAMES = {"Quad4": "Q4", "Tri3": "Tri3", "Hex8": "Hex8", "Tet4": "Tet4",
                 "Tet10": "Tet10"}
ELEMENT_DOMAINS = {"Quad4": "Q4 plate", "Hex8": "Hex8 block", "Tet4": "Tet4 block"}

#: sparlab_bench file stem per (element, solver), as scripts/run_scaling.sh
#: writes them.
SCALING_STEMS = {
    ("Quad4", "ldlt"): "runtime_scaling",
    ("Quad4", "amg"): "runtime_scaling_amg_cg",
    ("Hex8", "ldlt"): "runtime_scaling_3d",
    ("Hex8", "amg"): "runtime_scaling_3d_amg_cg",
    ("Hex8", "jacobi"): "runtime_scaling_3d_conjugate_gradient",
    ("Tet4", "ldlt"): "runtime_scaling_3d_tet",
    ("Tet4", "amg"): "runtime_scaling_3d_tet_amg_cg",
}


def load_solver_scaling(directory: str) -> Dict[Tuple[str, str], Tuple[pd.DataFrame, Dict]]:
    """Every solver-scaling table in `directory`, keyed by (element, solver)."""
    tables: Dict[Tuple[str, str], Tuple[pd.DataFrame, Dict]] = {}
    for key, stem in SCALING_STEMS.items():
        csv_path = os.path.join(directory, f"{stem}.csv")
        if not os.path.isfile(csv_path):
            continue
        table = load_csv(csv_path).sort_values("num_dofs").reset_index(drop=True)
        json_path = os.path.join(directory, f"{stem}.json")
        meta = load_json(json_path) if os.path.isfile(json_path) else {}
        tables[key] = (table, meta)
    return tables


def loglog_slope(x, y, tail: int = 3) -> float:
    """Least-squares slope of log(y) on log(x) over the last `tail` points."""
    x = np.asarray(x, dtype=float)[-tail:]
    y = np.asarray(y, dtype=float)[-tail:]
    keep = (x > 0) & (y > 0)
    if keep.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(x[keep]), np.log(y[keep]), 1)[0])


def _log_log_axes(ax) -> None:
    """Log scales on both axes, labelling only the decades on an axis that
    spans one or more: minor-tick labels crowd a small panel."""
    from matplotlib.ticker import NullFormatter

    ax.set_xscale("log")
    ax.set_yscale("log")
    for axis, (low, high) in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
        if low > 0 and high / low >= 10.0:
            axis.set_minor_formatter(NullFormatter())


def _one_two_five(axis) -> None:
    """Label a log axis at 1, 2 and 5 times each decade, as plain numbers: for
    counts (iterations, entries per DOF) that span one or two decades."""
    from matplotlib.ticker import LogLocator, NullFormatter, StrMethodFormatter

    axis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 5.0)))
    axis.set_major_formatter(StrMethodFormatter("{x:g}"))
    axis.set_minor_formatter(NullFormatter())


def _grid_with_legend_row(width: float, height: float, nrows: int, ncols: int):
    """A `nrows` x `ncols` panel grid plus a thin full-width row beneath it for a
    legend shared by every panel (colour means the same thing in all of them).

    A figure-level legend placed "outside" would compete with the footnote for
    the bottom margin; a row of its own is laid out like any other panel."""
    import matplotlib.pyplot as plt

    st.apply_style()
    fig = plt.figure(figsize=(width, height), layout="constrained")
    spec = fig.add_gridspec(nrows + 1, ncols, height_ratios=[1.0] * nrows + [0.16])
    grid = np.array([[fig.add_subplot(spec[r, c]) for c in range(ncols)]
                     for r in range(nrows)])
    legend_ax = fig.add_subplot(spec[nrows, :])
    legend_ax.axis("off")
    return fig, grid, legend_ax


def _solve_time(table: pd.DataFrame) -> np.ndarray:
    """Factorisation (or preconditioner setup) plus one solve, from scratch."""
    return (table["factorize[s]"] + table["solve[s]"]).to_numpy()


def _ratio_label(ratio: float) -> str:
    return f"{ratio:.1f}x" if ratio < 10.0 else f"{ratio:.0f}x"


def plot_solver_scaling(directory: str, path: str) -> str:
    """Time to solve from scratch and solver storage, per solver and element."""
    tables = load_solver_scaling(directory)
    if not tables:
        raise FileNotFoundError(
            f"no runtime_scaling*.csv in {directory}; run scripts/run_scaling.sh")
    elements = [e for e in ("Quad4", "Hex8", "Tet4") if any(k[0] == e for k in tables)]
    meta = next(iter(tables.values()))[1]

    fig, grid, legend_ax = _grid_with_legend_row(7.8, 6.4, 2, len(elements))
    for column, element in enumerate(elements):
        ax_time, ax_store = grid[0, column], grid[1, column]
        marker = ELEMENT_MARKERS[element]
        reference = None
        for solver in ("ldlt", "amg", "jacobi"):
            entry = tables.get((element, solver))
            if entry is None:
                continue
            table = entry[0]
            dofs = table["num_dofs"].to_numpy(dtype=float)
            total = _solve_time(table)
            color = st.series_color(SOLVER_SLOTS[solver])
            ax_time.plot(dofs, total, marker + "-", color=color, markersize=4.2,
                         linewidth=1.6, label=SOLVER_NAMES[solver])
            # The fitted slope as a direct label at the line end.
            ax_time.annotate(f"{loglog_slope(dofs, total):.2f}", (dofs[-1], total[-1]),
                             textcoords="offset points", xytext=(4, 0), fontsize=7.2,
                             color=st.INK_SECONDARY, va="center", ha="left")
            if solver != "jacobi":
                ax_store.plot(dofs, table["solver_nonzeros"].to_numpy(dtype=float) / dofs,
                              marker + "-", color=color, markersize=4.2, linewidth=1.6,
                              label=SOLVER_NAMES[solver])
            if reference is None or len(table) > len(reference):
                reference = table

        # Speed-up at the largest size both the direct and the multigrid solver ran.
        ldlt, amg = tables.get((element, "ldlt")), tables.get((element, "amg"))
        if ldlt is not None and amg is not None:
            both = pd.merge(ldlt[0], amg[0], on="num_dofs", suffixes=("_d", "_m"))
            if not both.empty:
                last = both.iloc[-1]
                t_direct = float(last["factorize[s]_d"] + last["solve[s]_d"])
                t_multigrid = float(last["factorize[s]_m"] + last["solve[s]_m"])
                ax_time.plot([last["num_dofs"]] * 2, [t_multigrid, t_direct], "-",
                             color=st.INK_MUTED, linewidth=0.9)
                ax_time.annotate(_ratio_label(t_direct / t_multigrid),
                                 (last["num_dofs"], np.sqrt(t_direct * t_multigrid)),
                                 textcoords="offset points", xytext=(-4, 0), fontsize=7.6,
                                 color=st.INK_PRIMARY, va="center", ha="right")
        if reference is not None:
            dofs = reference["num_dofs"].to_numpy(dtype=float)
            ax_store.plot(dofs, reference["stiffness_nonzeros"].to_numpy(dtype=float) / dofs,
                          "--", color=st.INK_MUTED, linewidth=1.1,
                          label="the stiffness matrix itself")

        for ax in (ax_time, ax_store):
            _log_log_axes(ax)
        _one_two_five(ax_store.yaxis)
        ax_time.set_title(ELEMENT_DOMAINS[element], loc="left", fontsize=10)
        ax_store.set_xlabel("degrees of freedom")
        if column == 0:
            ax_time.set_ylabel("factorise or set up, then solve [s]")
            ax_store.set_ylabel("stored entries per DOF [-]")

    # One legend per row for the whole figure: the colour is the solver in every
    # panel, and the marker (the element) is named by the panel title.
    from matplotlib.lines import Line2D

    present = [s for s in ("ldlt", "amg", "jacobi") if any(k[1] == s for k in tables)]
    handles = [Line2D([], [], color=st.series_color(SOLVER_SLOTS[s]), linewidth=1.8,
                      label=SOLVER_NAMES[s]) for s in present]
    handles.append(Line2D([], [], color=st.INK_MUTED, linewidth=1.1, linestyle="--",
                          label="entries of K itself (lower row)"))
    legend_ax.legend(handles=handles, loc="center", ncol=2, fontsize=7.8, frameon=False)

    threads = meta.get("openmp_threads")
    st.figure_title(
        fig, "Linear solvers: cost of one solve from scratch, and storage",
        f"{meta.get('compiler', 'unknown compiler')}, {meta.get('build_type', '')}"
        + (f", {threads} OpenMP threads" if threads else "")
        + "; CG to a relative residual of 1e-10; minimum over repeats. Numbers at "
        "the line ends are least-squares slopes of log(time) on log(DOFs) over the "
        "largest three sizes. The vertical bar is labelled with the Cholesky time "
        "over the multigrid time at the largest size both solvers ran.",
    )
    st.annotate_note(
        fig,
        "Upper row: sparse Cholesky is AMD-ordered factorisation plus one "
        "back-substitution; the CG solvers are preconditioner setup plus the "
        "iterations for one load case from a zero initial guess. Lower row: "
        "entries the solver stores beyond K itself - the factor L and D for "
        "Cholesky; the coarse operators, prolongators, restrictions and dense "
        "coarsest factor for multigrid. Jacobi stores only a diagonal. The "
        "benchmark problem is a cantilever with a tip load (sparlab_bench).",
    )
    return st.save_figure(fig, path)


def plot_solver_iterations(directory: str, path: str) -> str:
    """CG iteration counts against problem size for both preconditioners."""
    tables = load_solver_scaling(directory)
    series = [(k, v) for k, v in tables.items() if k[1] in ("amg", "jacobi")]
    if not series:
        raise FileNotFoundError(
            f"no iterative runtime_scaling*.csv in {directory}; run scripts/run_scaling.sh")
    order = {("Quad4", "amg"): 0, ("Hex8", "amg"): 1, ("Tet4", "amg"): 2,
             ("Hex8", "jacobi"): 3}
    series.sort(key=lambda item: order.get(item[0], 9))

    fig, ax = st.figure(7.2, 3.9)
    single_level = False
    for (element, solver), (table, _meta) in series:
        color = st.series_color(SOLVER_SLOTS[solver])
        marker = ELEMENT_MARKERS[element]
        dofs = table["num_dofs"].to_numpy(dtype=float)
        iterations = table["iterations"].to_numpy(dtype=float)
        # A multigrid "hierarchy" of one level is the direct coarse solve: one
        # iteration by construction, drawn hollow and left out of the fit.
        multi = (table["levels"].to_numpy() >= 2) if solver == "amg" else np.ones(len(table), bool)
        slope = loglog_slope(dofs[multi], iterations[multi], tail=len(dofs))
        label = f"{SOLVER_NAMES[solver]}, {ELEMENT_NAMES[element]}"
        if np.isfinite(slope):
            label += f" (slope {slope:.2f})"
        ax.plot(dofs[multi], iterations[multi], marker + "-", color=color, markersize=5,
                linewidth=1.6, label=label)
        if not multi.all():
            single_level = True
            ax.plot(dofs[~multi], iterations[~multi], marker, color=color, markersize=5,
                    markerfacecolor="none", markeredgewidth=1.2)
    _log_log_axes(ax)
    _one_two_five(ax.yaxis)
    ax.set_xlabel("degrees of freedom")
    ax.set_ylabel("CG iterations to 1e-10 [-]")
    st.title(
        ax, "Conjugate-gradient iterations under refinement",
        "slopes are least-squares fits of log(iterations) on log(DOFs) over every "
        "multi-level size; a scalable preconditioner keeps the count flat",
    )
    st.legend(ax, loc="lower right", fontsize=7.8)
    note = ("Colour is the preconditioner, marker the element. Jacobi CG needs "
            "O(sqrt(condition number)) iterations, which grows as 1/h: n^(1/3) in 3-D "
            "and n^(1/2) in 2-D.")
    if single_level:
        note += (" Hollow markers are sizes below the coarse-grid size, where the "
                 "multigrid solver is a direct solve and converges in one iteration.")
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: the 3-D (Hex8) studies
# ---------------------------------------------------------------------------
def plot_sensitivity_check_3d(directory: str, path: str) -> str:
    """Error of the Hex8 topology gradient vs the central-difference step."""
    steps = load_csv(os.path.join(directory, "sensitivity_steps_3d.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get("sensitivity_3d", {})
    tolerance = summary.get("tolerance", 1e-5)

    fig, ax = st.figure(7.2, 3.6)
    st.require_scatter_series(3, "error measures")
    for slot, (column, label) in enumerate(
        [
            ("max_relative_error[-]", "max relative error over elements"),
            ("rms_relative_error[-]", "RMS relative error over elements"),
            ("directional_relative_error[-]", "relative error of the directional derivative"),
        ]
    ):
        ax.plot(steps["step[-]"], steps[column], "o-", color=st.series_color(slot),
                label=label)
    ax.axhline(tolerance, color=st.INK_MUTED, linewidth=1.0, linestyle="--")
    ax.text(steps["step[-]"].min(), tolerance, f" pass tolerance {tolerance:g}",
            fontsize=7.8, color=st.INK_SECONDARY, va="bottom")
    _log_log_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("central-difference step h on the design variable [-]")
    ax.set_ylabel("relative error [-]")
    st.title(
        ax, "Hex8 topology sensitivity: analytical gradient vs central differences",
        f"{summary.get('mesh', 'Hex8 mesh')}, {summary.get('filter', 'density filter')}; "
        f"{summary.get('tested_elements', '?')} elements tested, "
        f"{summary.get('excluded_elements', 0)} excluded at active bounds; best step "
        f"{summary.get('best_step', float('nan')):g} gives "
        f"{summary.get('best_max_relative_error', float('nan')):.2e}",
    )
    st.legend(ax, loc="upper left")
    st.annotate_note(
        fig,
        "The same adjoint gradient as in 2-D, evaluated with the trilinear "
        "element; the error floor is where round-off in the difference of two "
        "compliances overtakes the truncation error of the central difference.",
    )
    return st.save_figure(fig, path)


def plot_mesh_convergence_3d(directory: str, path: str) -> str:
    """Hex8 cantilever tip deflection and discretisation error vs element size."""
    # Coarse to fine, so the last row is the finest mesh the self-convergence
    # error is measured against.
    table = load_csv(os.path.join(directory, "mesh_convergence_3d.csv")).sort_values(
        "h[m]", ascending=False)
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("mesh_convergence_3d", {})
    order = block.get("observed_convergence_order_tip_deflection")

    fig, axes = st.figure(7.4, 6.4, nrows=2, ncols=1)

    ax = axes[0]
    ax.plot(table["h[m]"], np.abs(table["tip_mean[m]"]), "o-", color=st.series_color(0),
            label="FEM, Hex8, full 2x2x2 integration")
    ax.axhline(abs(float(table["timoshenko[m]"].iloc[0])), color=st.series_color(0),
               linewidth=0.9, linestyle=":", label="Timoshenko (bending + shear)")
    ax.axhline(abs(float(table["euler_bernoulli[m]"].iloc[0])), color=st.INK_MUTED,
               linewidth=1.0, linestyle="--", label="Euler-Bernoulli (bending only)")
    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("|tip deflection| [m]")
    st.title(
        ax, "Solid cantilever tip deflection vs mesh size",
        f"{block.get('geometry', '')}; deflection is the mean u_y over the tip face",
    )
    st.legend(ax, loc="lower left")

    ax = axes[1]
    finest = float(table["tip_mean[m]"].iloc[-1])
    error = np.abs(table["tip_mean[m]"] - finest) / abs(finest)
    mask = error > 0
    label = "self-convergence error"
    if order is not None:
        label += f" (observed order {order:.2f})"
    ax.plot(table["h[m]"][mask], error[mask], "o-", color=st.series_color(0), label=label)
    h_ref = np.array([float(table["h[m]"].min()) * 2.0, float(table["h[m]"].max())])
    ax.plot(h_ref, 1e-3 * (h_ref / h_ref[0]) ** 2, "-", color=st.INK_MUTED,
            linewidth=1.0, label=r"reference slope $h^2$")
    _log_log_axes(ax)
    ax.invert_xaxis()
    ax.set_xlabel("element size h [m]")
    ax.set_ylabel("relative error vs the finest mesh [-]")
    st.title(
        ax, "Self-convergence of the discretisation error",
        f"finest mesh within {block.get('finest_mesh_relative_error_vs_timoshenko', float('nan')):.2e} "
        "of Timoshenko",
    )
    st.legend(ax, loc="lower left")
    st.annotate_note(
        fig,
        block.get("note", "") + " Comparing against the finest mesh is verification; "
        "comparing against beam theory is validation, where a finite gap is expected.",
    )
    return st.save_figure(fig, path)


def plot_modal_convergence_3d(directory: str, path: str) -> str:
    """Hex8 cantilever frequencies for the weak and strong bending axes."""
    table = load_csv(os.path.join(directory, "modal_convergence_3d.csv")).sort_values("num_dofs")
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("modal_validation_3d", {})
    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.5)

    ax = axes[0]
    series = [
        ("f1_weak_fem[Hz]", "f1_weak_theory[Hz]", "o", "first bending mode, weak axis"),
        ("f1_strong_fem[Hz]", "f1_strong_theory[Hz]", "s", "first bending mode, strong axis"),
    ]
    for slot, (column, reference, marker, label) in enumerate(series):
        ax.plot(table["num_dofs"], table[column], marker + "-", color=st.series_color(slot),
                label=f"{label}, FEM")
        ax.axhline(float(table[reference].iloc[0]), color=st.series_color(slot),
                   linewidth=0.9, linestyle="--", label=f"{label}, Euler-Bernoulli")
    ax.set_xscale("log")
    ax.set_ylabel("frequency [Hz]")
    st.title(
        ax, "Solid cantilever natural frequencies vs mesh size",
        f"{block.get('geometry', '')}; consistent mass matrix; mass conserved to "
        f"{block.get('max_mass_relative_error', float('nan')):.1e} relative",
    )
    st.legend(ax, loc="center right", ncol=2)

    ax = axes[1]
    for slot, (column, marker, label) in enumerate([
        ("f1_weak_rel_error[-]", "o", "weak axis vs Euler-Bernoulli"),
        ("f1_strong_rel_error[-]", "s", "strong axis vs Euler-Bernoulli"),
    ]):
        ax.plot(table["num_dofs"], table[column], marker + "-", color=st.series_color(slot),
                label=label)
    _log_log_axes(ax)
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("relative difference [-]")
    st.title(ax, "difference from beam theory", block.get("note", ""))
    st.legend(ax, loc="upper right")
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Verification: linear simplices, multigrid, sensitivities through projection
# ---------------------------------------------------------------------------
def plot_mesh_convergence_simplex(directory: str, path: str) -> str:
    """Tri3 / Tet4 cantilever error against the Q4 / Hex8 results on the same beams."""
    simplex = load_csv(os.path.join(directory, "mesh_convergence_simplex.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    block = summary.get("mesh_convergence_simplex", {})

    panels = []
    q4_path = os.path.join(directory, "mesh_convergence.csv")
    q4 = load_csv(q4_path) if os.path.isfile(q4_path) else None
    if q4 is not None:
        q4 = q4[np.isclose(q4["poisson[-]"], 0.3)]
    q4_order = (summary.get("mesh_convergence", {})
                .get("observed_convergence_order_tip_deflection", {}).get("0.3"))
    panels.append(("Tri3", q4, "Q4", q4_order,
                   "plane-stress cantilever L = 1 m, h = 0.1 m, t = 10 mm, "
                   "E = 70 GPa, nu = 0.3, 1 kN tip load"))
    hex_path = os.path.join(directory, "mesh_convergence_3d.csv")
    hexa = load_csv(hex_path) if os.path.isfile(hex_path) else None
    hex_order = summary.get("mesh_convergence_3d", {}).get(
        "observed_convergence_order_tip_deflection")
    panels.append(("Tet4", hexa, "Hex8", hex_order,
                   "solid cantilever 1 x 0.1 x 0.05 m, E = 70 GPa, nu = 0, "
                   "1 kN tip load"))

    fig, axes = st.figure(7.4, 6.4, nrows=2, ncols=1)
    for ax, (element, reference, ref_name, ref_order, geometry) in zip(axes, panels):
        sub = simplex[simplex["element"] == element].sort_values("num_dofs")
        if sub.empty:
            ax.axis("off")
            continue
        order = block.get(element, {}).get("observed_convergence_order_tip_deflection")
        label = f"{element} (linear simplex)"
        if order is not None:
            label += f", observed order {order:.2f}"
        ax.plot(sub["num_dofs"], sub["rel_error_timoshenko[-]"], "o-",
                color=st.series_color(0), label=label)
        if reference is not None and not reference.empty:
            reference = reference.sort_values("num_dofs")
            label = f"{ref_name} on the same beam"
            if ref_order is not None:
                label += f", observed order {ref_order:.2f}"
            ax.plot(reference["num_dofs"], reference["rel_error_timoshenko[-]"], "s-",
                    color=st.series_color(1), label=label)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("number of degrees of freedom")
        ax.set_ylabel("relative error vs Timoshenko [-]")
        finest = block.get(element, {}).get("finest_mesh_relative_error_vs_timoshenko")
        st.title(
            ax, f"{element}: tip deflection error vs mesh size",
            geometry + (f"; finest {element} mesh within {100 * finest:.2f} % of "
                        "Timoshenko" if finest is not None else ""),
        )
        st.legend(ax, loc="upper right")
    st.annotate_note(
        fig,
        "The simplex meshes split each structured cell (2 triangles per quadrilateral, "
        "6 tetrahedra per hexahedron), so both curves in a panel describe the same "
        "beam. Observed orders are in element size h, from self-convergence over "
        "the three finest meshes. Constant-strain simplices lock in bending and "
        "approach the beam value from below, more slowly per DOF than the bilinear "
        "and trilinear elements. The error is measured against beam theory, so it "
        "also contains the modelling gap between a beam and plane or solid "
        "elasticity, which is where the Q4 curve levels off.",
    )
    return st.save_figure(fig, path)


def plot_multigrid_study(directory: str, path: str) -> str:
    """Multigrid CG against Cholesky and Jacobi CG as the mesh is refined."""
    table = load_csv(os.path.join(directory, "multigrid_scaling.csv"))
    block = load_json(os.path.join(directory, "summary.json")).get("multigrid", {})
    elements = [e for e in ("Hex8", "Tet4") if (table["element"] == e).any()]
    if not elements:
        raise ValueError("multigrid_scaling.csv has no Hex8 or Tet4 rows")

    fig, grid, legend_ax = _grid_with_legend_row(7.6, 6.1, 2, len(elements))
    single_level = False
    for column, element in enumerate(elements):
        sub = table[table["element"] == element].sort_values("num_dofs")
        dofs = sub["num_dofs"].to_numpy(dtype=float)
        multi = sub["amg_levels"].to_numpy() >= 2
        marker = ELEMENT_MARKERS[element]
        amg, jacobi = (st.series_color(SOLVER_SLOTS[s]) for s in ("amg", "jacobi"))

        ax = grid[0, column]
        ax.plot(dofs, sub["jacobi_iterations"], marker + "-", color=jacobi, markersize=4.5,
                label=SOLVER_NAMES["jacobi"])
        ax.plot(dofs[multi], sub["amg_iterations"].to_numpy()[multi], marker + "-",
                color=amg, markersize=4.5, label=SOLVER_NAMES["amg"])
        if not multi.all():
            single_level = True
            ax.plot(dofs[~multi], sub["amg_iterations"].to_numpy()[~multi], marker,
                    color=amg, markersize=4.5, markerfacecolor="none", markeredgewidth=1.2)
        for x, y, levels in zip(dofs[multi], sub["amg_iterations"].to_numpy()[multi],
                                sub["amg_levels"].to_numpy()[multi]):
            ax.annotate(f"{int(levels)} levels", (x, y), textcoords="offset points",
                        xytext=(0, -11), ha="center", fontsize=6.8, color=st.INK_SECONDARY)
        _log_log_axes(ax)
        _one_two_five(ax.yaxis)
        ax.set_title(f"{element}: CG iterations to 1e-10", loc="left", fontsize=10)
        if column == 0:
            ax.set_ylabel("iterations [-]")

        ax = grid[1, column]
        for solver, key in (("ldlt", "ldlt_seconds[s]"), ("amg", "amg_seconds[s]"),
                            ("jacobi", "jacobi_seconds[s]")):
            ax.plot(dofs, sub[key], marker + "-", markersize=4.5,
                    color=st.series_color(SOLVER_SLOTS[solver]), label=SOLVER_NAMES[solver])
        _log_log_axes(ax)
        ax.set_xlabel("degrees of freedom")
        ax.set_title(f"{element}: wall-clock time of one solve", loc="left", fontsize=10)
        if column == 0:
            ax.set_ylabel("seconds [s]")
    from matplotlib.lines import Line2D

    legend_ax.legend(
        handles=[Line2D([], [], color=st.series_color(SOLVER_SLOTS[s]), linewidth=1.8,
                        label=SOLVER_NAMES[s]) for s in ("ldlt", "amg", "jacobi")],
        loc="center", ncol=3, fontsize=7.8, frameon=False)

    subtitle = ("cantilever block 2 x 1 x 0.5 m (nx x nx/2 x nx/4 cells), E = 70 GPa, "
                "nu = 0.3, tip load")
    difference = block.get("max_relative_difference_vs_ldlt")
    if difference is None:
        difference = float(table["relative_difference_vs_ldlt[-]"].max())
    subtitle += ("; largest relative displacement difference from the Cholesky "
                 f"solution {difference:.1e}")
    growth = block.get("max_iteration_growth_coarsest_to_finest")
    if growth is not None:
        subtitle += ("; iteration growth from the coarsest to the finest multi-level "
                     f"mesh {growth:.2f}x")
    st.figure_title(
        fig, "Multigrid CG: agreement with Cholesky and iterations under refinement",
        subtitle,
    )
    note = ("Times include the preconditioner setup or the factorisation. "
            "The number under each multigrid point is the depth of its hierarchy.")
    if single_level:
        note += (" Hollow markers are meshes below the coarse-grid size, where the "
                 "multigrid solver is a direct solve and converges in one iteration.")
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


def plot_sensitivity_projection(directory: str, path: str) -> str:
    """Gradient error through the Heaviside projection vs the finite-difference step."""
    table = load_csv(os.path.join(directory, "sensitivity_projection.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    tolerance = next((o["tolerance"] for o in summary.get("outcomes", [])
                      if "Heaviside" in o.get("study", "")), None)
    elements = [e for e in ("Quad4", "Tet4") if (table["element"] == e).any()]
    betas = sorted(table["beta[-]"].unique())
    # beta is ordinal, so a sequential ramp; its light end stops short of yellow
    # so the thin lines keep their contrast on the page.
    colors = st.sequence_colors(len(betas), lo=0.08, hi=0.68)

    fig, axes, legend_ax = _grid_with_legend_row(7.4, 6.9, len(elements), 1)
    for row, element in enumerate(elements):
        ax = axes[row, 0]
        for color, beta in zip(colors, betas):
            sub = table[(table["element"] == element)
                        & np.isclose(table["beta[-]"], beta)].sort_values("step[-]")
            ax.plot(sub["step[-]"], sub["max_scaled_error[-]"], "o-", color=color)
            ax.plot(sub["step[-]"], sub["directional_relative_error[-]"], "s--",
                    color=color, markersize=4, linewidth=1.3)
        if tolerance:
            ax.axhline(tolerance, color=st.INK_MUTED, linewidth=1.0, linestyle=":")
            ax.text(float(table["step[-]"].min()), tolerance, f"pass tolerance {tolerance:g}",
                    fontsize=7.6, color=st.INK_SECONDARY, va="bottom", ha="right")
        _log_log_axes(ax)
        ax.invert_xaxis()
        ax.set_xlabel("central-difference step h on the design variable [-]")
        ax.set_ylabel("relative error [-]")
        tested = int(table[table["element"] == element]["num_tested"].iloc[0])
        st.title(
            ax, f"{ELEMENT_NAMES[element]}: compliance gradient through filter and projection",
            f"{tested} design variables tested at eta = 0.5 on a smoothly varying "
            "design; the chain rule runs through the Heaviside projection and the "
            "density filter",
        )
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color=color, linewidth=1.9, label=f"beta = {beta:g}")
               for color, beta in zip(colors, betas)]
    handles += [
        Line2D([], [], color=st.INK_SECONDARY, marker="o", linewidth=1.9,
               label="largest entry error"),
        Line2D([], [], color=st.INK_SECONDARY, marker="s", markersize=4, linewidth=1.3,
               linestyle="--", label="directional derivative"),
    ]
    legend_ax.legend(handles=handles, loc="center", ncol=len(handles), fontsize=7.6,
                     frameon=False, handlelength=2.2, columnspacing=1.2)
    st.annotate_note(
        fig,
        "Each entry is judged against max(|analytical|, |finite difference|, "
        "1e-3 ||gradient||_inf): at large beta the projection's derivative "
        "vanishes away from the threshold, and there the comparison is with the "
        "round-off floor of a central difference rather than with a relative "
        "error of a number near zero. The study passes on the best step of each "
        "beta, taking the worse of the two measures.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Cross-validation against independent codes
# ---------------------------------------------------------------------------
def plot_cross_validation(directory: str, path: str) -> str:
    """Largest relative nodal-displacement difference per code and load case."""
    summary = load_json(os.path.join(directory, "summary.json"))
    codes = summary.get("codes", {})
    tolerances = summary.get("tolerances", {})
    rows = []
    for case in summary.get("cases", []):
        for load_case in case.get("load_cases", []):
            rows.append((case["case"], case.get("element_type", ""),
                         load_case["load_case"], load_case["codes"]))
    if not rows:
        raise ValueError("cross-validation summary lists no comparisons")

    # The slot follows the code, not its position among the codes present, so
    # scikit-fem keeps its colour on a runner where CalculiX is not installed.
    code_names = ["scikit-fem", "calculix"]
    st.require_scatter_series(len(code_names), "codes")
    present = [code for code in code_names if any(code in r[3] for r in rows)]
    # The legend gets a row of its own beneath the panel: beside it, it took a
    # third of the width from a log axis spanning ten decades.
    fig, grid, legend_ax = _grid_with_legend_row(7.4, 0.36 * len(rows) + 3.0, 1, 1)
    ax = grid[0, 0]
    y = np.arange(len(rows))[::-1]
    floors = set()
    informational = 0
    for slot, code in enumerate(code_names):
        if code not in present:
            continue
        xs, ys, info_x, info_y = [], [], [], []
        for position, (_c, _e, _l, results) in zip(y, rows):
            entry = results.get(code)
            if entry is None:
                continue
            # passed is None for a comparison between two different
            # idealisations (a plane element CalculiX expands through the
            # thickness, at nu != 0): recorded, not judged.
            if entry.get("passed") is None:
                info_x.append(entry["max_rel_diff"])
                info_y.append(position)
            else:
                xs.append(entry["max_rel_diff"])
                ys.append(position)
            floor = entry.get("frd_rounding_floor_rel")
            if floor:
                floors.add(float(floor))
        version = codes.get(code, {}).get("version", "")
        ax.plot(xs, ys, "o", color=st.series_color(slot), markersize=7,
                label=f"{code} {version}".strip())
        if info_x:
            informational += len(info_x)
            ax.plot(info_x, info_y, "o", color=st.series_color(slot), markersize=7,
                    markerfacecolor="none", markeredgewidth=1.5,
                    label=f"{code}: different idealisation, not judged")
    tol_lines = [
        ("skfem", "scikit-fem tolerance", 0, "--", "scikit-fem"),
        ("calculix_solid", "CalculiX tolerance", 1, "--", "calculix"),
    ]
    for key, label, slot, style, code in tol_lines:
        value = tolerances.get(key)
        if value and code in present:
            ax.axvline(value, color=st.series_color(slot), linewidth=1.0, linestyle=style,
                       label=f"{label} {value:g}")
    for floor in sorted(floors):
        ax.axvline(floor, color=st.INK_MUTED, linewidth=1.0, linestyle=":",
                   label=f".frd six-digit rounding floor {floor:g}")
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} ({e}), '{l}'" for c, e, l, _r in rows], fontsize=8.0)
    ax.set_xlabel("max |u_SparLab - u_reference| / max |u_reference| over all nodes [-]")
    ax.set_ylim(-0.7, len(rows) - 0.3)
    judged = [entry["passed"] for *_r, results in rows for entry in results.values()
              if entry.get("passed") is not None]
    verdict = ("every judged comparison within its tolerance" if all(judged)
               else "a comparison FAILED")
    if informational:
        verdict += f"; {informational} informational"
    st.figure_title(
        fig, "Cross-validation: nodal displacements vs independent codes",
        f"{len(rows)} load cases, {len(present)} "
        f"code{'s' if len(present) != 1 else ''}; {verdict}",
    )
    handles, labels = ax.get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="center", ncol=2, fontsize=7.8, frameon=False,
                     columnspacing=1.6)
    st.annotate_note(
        fig,
        "Same mesh, material, supports and nodal loads in every code. scikit-fem "
        "uses the same elements (bilinear, trilinear and linear simplices), so its "
        "differences are solver round-off. CalculiX C3D8 and C3D4 are the same "
        "elements as SparLab's Hex8 and Tet4. Its plane elements (CPS4, CPS3) are "
        "expanded into a layer of solid elements, which matches plane stress only "
        "for nu = 0; plane comparisons at nu != 0 are drawn hollow and not judged. "
        "CalculiX results are read from the .frd file, which carries six "
        "significant digits, so differences below the dotted floor are its output "
        "rounding.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Heaviside projection: grey level and the SIMP-vs-structure gap
# ---------------------------------------------------------------------------
#: (label, run without projection, run with projection). Either may be None.
PROJECTION_RUNS = [
    ("MBB beam (2-D, Q4)", "mbb_beam", "mbb_beam_projected"),
    ("bracket (3-D, Hex8)", "bracket_3d", "bracket_3d_projected"),
    ("lug bracket (2-D, Tri3, Gmsh)", None, "lug_bracket_2d"),
    ("engine mount (3-D, Tet4, Gmsh)", None, "engine_mount_3d"),
    ("large bracket (3-D, Hex8, 356k DOFs)", None, "bracket_3d_large"),
]


def collect_projection(results_dir: str) -> pd.DataFrame:
    """One row per run in PROJECTION_RUNS that exists under `results_dir`."""
    rows = []
    for label, plain, projected in PROJECTION_RUNS:
        for variant, name in (("without projection", plain), ("with projection", projected)):
            if name is None:
                continue
            path = os.path.join(results_dir, name, "summary.json")
            if not os.path.isfile(path):
                continue
            doc = load_json(path)
            result = doc.get("optimization_result", {})
            solid = doc.get("interpreted_solid_analysis") or {}
            projection = result.get("projection") or {}
            rows.append({
                "problem": label,
                "run": name,
                "variant": variant,
                "final_beta": projection.get("final_beta"),
                "iterations": result.get("iterations"),
                "stop_reason": result.get("stop_reason"),
                "grey_level": result.get("grey_level"),
                "filtered_grey_level": projection.get("filtered_grey_level"),
                "compliance_J": result.get("compliance_J"),
                "interpreted_compliance_J": solid.get("weighted_compliance_J"),
                "compliance_vs_simp_ratio": solid.get("compliance_vs_simp_ratio"),
                "seconds": result.get("total_seconds"),
            })
    if not rows:
        raise FileNotFoundError(
            f"none of the projection runs exist under {results_dir}; run "
            "scripts/run_all_benchmarks.sh")
    return pd.DataFrame(rows)


def plot_projection_summary(results_dir: str, path: str) -> str:
    """Grey level and thresholded-vs-SIMP compliance, with and without projection."""
    table = collect_projection(results_dir)
    problems = [label for label, *_ in PROJECTION_RUNS if (table["problem"] == label).any()]
    fig, axes = st.figure(7.6, 0.42 * len(problems) + 3.0, nrows=1, ncols=2, sharey=True)
    y = {label: position for position, label in enumerate(problems[::-1])}
    styles = {"without projection": (0, "o"), "with projection": (1, "D")}
    for ax, column, xlabel in (
        (axes[0], "grey_level", "grey level [-]"),
        (axes[1], "compliance_vs_simp_ratio", "thresholded / SIMP compliance [-]"),
    ):
        for label in problems:
            sub = table[table["problem"] == label]
            if len(sub) == 2 and sub[column].notna().all():
                ax.plot(sub[column], [y[label]] * 2, "-", color=st.GRID, linewidth=3.0,
                        zorder=1)
        for variant, (slot, marker) in styles.items():
            sub = table[table["variant"] == variant]
            ax.plot(sub[column], [y[p] for p in sub["problem"]], marker,
                    color=st.series_color(slot), markersize=7, zorder=3,
                    label=variant)
            for _, row in sub.iterrows():
                if row[column] is None or row[column] != row[column]:
                    continue
                ax.annotate(f"{row[column]:.3f}", (row[column], y[row["problem"]]),
                            textcoords="offset points", xytext=(0, 7), ha="center",
                            fontsize=7.0, color=st.INK_SECONDARY)
        ax.set_xlabel(xlabel, fontsize=8.6)
        ax.set_ylim(-0.7, len(problems) - 0.2)
    axes[1].axvline(1.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
    ratios = table["compliance_vs_simp_ratio"].dropna()
    if not ratios.empty:
        axes[1].set_xlim(min(0.9, float(ratios.min()) - 0.05),
                         max(1.05, float(ratios.max()) + 0.03))
    greys = table["grey_level"].dropna()
    axes[0].set_xlim(0.0, max(0.05, float(greys.max()) * 1.15) if not greys.empty else 1.0)
    axes[0].set_yticks(list(y.values()))
    axes[0].set_yticklabels(list(y.keys()), fontsize=8.2)
    st.legend(axes[0], loc="lower right", fontsize=7.8)
    st.figure_title(
        fig, "Heaviside projection: how binary the design is, and what that buys",
        "grey level = 4 mean(rho (1 - rho)) over all elements; the ratio "
        "re-solves the design thresholded at 0.5 as solid material, so 1.0 means "
        "the optimiser's objective is the structure's real compliance",
    )
    st.annotate_note(
        fig,
        "The grey rails join the two runs of one problem; the runs differ only in "
        "the projection block and the settings listed in each deck's header "
        "comment. Problems with one point were run with projection only. Values "
        "are read from each run's summary.json.",
    )
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Design study
# ---------------------------------------------------------------------------
def collect_study(study_dir: str) -> pd.DataFrame:
    """Flatten every `summary.json` under a study directory into one table."""
    rows: List[Dict] = []
    pattern = os.path.join(study_dir, "*", "*", "summary.json")
    for path in sorted(glob.glob(pattern)):
        parts = path.split(os.sep)
        arm, tag = parts[-3], parts[-2]
        doc = load_json(path)
        setup = doc.get("optimization_setup", {})
        result = doc.get("optimization_result", {})
        mesh = doc.get("mesh", {})
        material = doc.get("material", {})
        comparison = doc.get("mass_stiffness_comparison", {})
        interp = doc.get("solid_interpretation", {})
        row = {
            "arm": arm,
            "tag": tag,
            "volume_fraction_target": setup.get("volume_fraction_target"),
            "volume_fraction": result.get("volume_fraction"),
            "compliance_J": result.get("compliance_J"),
            "mass_kg": result.get("mass_kg"),
            "grey_level": result.get("grey_level"),
            "iterations": result.get("iterations"),
            "converged": result.get("converged"),
            "stop_reason": result.get("stop_reason"),
            "seconds": result.get("total_seconds"),
            "num_elements": mesh.get("num_elements"),
            "num_dofs": mesh.get("num_dofs"),
            "penalty": setup.get("simp_penalty"),
            "filter_radius_m": setup.get("filter_radius_m"),
            "filter_support": setup.get("filter_average_support"),
            "youngs_modulus_Pa": material.get("youngs_modulus_Pa"),
            "full_solid_compliance_J": comparison.get("full_solid_compliance_J"),
            "equal_mass_plate_compliance_J":
                comparison.get("equal_mass_uniform_plate_compliance_J"),
            "stiffness_gain_over_equal_mass_plate":
                comparison.get("stiffness_gain_over_equal_mass_plate"),
            "retained_elements": interp.get("elements_retained"),
            "connected_groups": interp.get("connected_groups_above_threshold"),
            "islands_discarded_m3": interp.get("volume_discarded_as_islands_m3"),
        }
        # What the thresholded structure actually delivers, which is not the
        # objective the optimiser minimised.
        solid = doc.get("interpreted_solid_analysis") or {}
        row["interpreted_compliance_J"] = solid.get("weighted_compliance_J")
        row["compliance_vs_simp_ratio"] = solid.get("compliance_vs_simp_ratio")
        row["interpreted_mass_kg"] = solid.get("mass_kg")
        row["interpreted_max_von_mises_Pa"] = solid.get("max_von_mises_Pa")
        row["interpreted_analysis_failed"] = bool(solid.get("analysis_failed", False))

        modal = doc.get("modal_optimised_topology")
        if modal and modal.get("frequencies_hz"):
            row["f1_topology_Hz"] = modal["frequencies_hz"][0]
            row["topology_mass_kg"] = modal.get("total_mass_kg")
        solid_modal = doc.get("modal_initial_solid")
        if solid_modal and solid_modal.get("frequencies_hz"):
            row["f1_solid_Hz"] = solid_modal["frequencies_hz"][0]
        # Per-load-case columns. The objective is the weight-normalised mean of
        # these, so both the weight and the case compliance are carried through:
        # compliance_J == sum_l weight_norm_l * compliance_<case_l>_J exactly.
        names = setup.get("load_case_names") or []
        for load_case, value in zip(names, result.get("load_case_compliance_J") or []):
            row[f"compliance_{load_case}_J"] = value
        for load_case, value in zip(names, setup.get("load_case_weights") or []):
            row[f"weight_{load_case}"] = value
        rows.append(row)

    if not rows:
        raise FileNotFoundError(
            f"no summary.json found under {study_dir}/*/*/; run "
            "scripts/run_aerospace_study.sh first"
        )
    return pd.DataFrame(rows)


def _numeric_tag(tag: str) -> float:
    """Extract the varied quantity from a study tag like `vf_0.35` or `p_3.0`.

    The *last* numeric token wins, not the first: a load-weighting tag is
    `w_<down>_<reversal>_<lateral>`, and the quantity the arm varies is the
    last one. A mesh tag's `AxB` token resolves to the element count, which is
    the right ordering key for `n_120x80`.

    Prefer a recorded column where one exists - `order_by` below - and keep this
    for the arms whose varied quantity is not a column of its own.
    """
    value = float("nan")
    for token in tag.replace("-", "_").split("_"):
        try:
            value = float(token)
            continue
        except ValueError:
            pass
        if "x" in token:
            parts = token.split("x")
            try:
                if len(parts) == 2:
                    value = float(parts[0]) * float(parts[1])
            except ValueError:
                pass
    return value


def plot_mass_stiffness_pareto(table: pd.DataFrame, path: str) -> str:
    """Compliance and first frequency against mass for the volume sweep."""
    arm = table[table["arm"] == "volume_fraction"].copy()
    if arm.empty:
        raise ValueError("the study table has no volume_fraction arm")
    arm["vf"] = arm["tag"].map(_numeric_tag)
    arm = arm.sort_values("mass_kg")

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.7)

    ax = axes[0]
    st.require_scatter_series(3, "mass-stiffness curves")
    ax.plot(arm["mass_kg"], arm["compliance_J"], "o-", color=st.series_color(0),
            label="SIMP density field (the objective minimised)")
    if "interpreted_compliance_J" in arm and arm["interpreted_compliance_J"].notna().any():
        # The structure a threshold actually produces. It is stiffer than the
        # objective here, because thresholding promotes the filter's grey
        # boundary band to solid material.
        ax.plot(arm["mass_kg"], arm["interpreted_compliance_J"], "^-",
                color=st.series_color(2),
                label="the same design thresholded at 0.5 and re-solved")
    ax.plot(arm["mass_kg"], arm["equal_mass_plate_compliance_J"], "s--",
            color=st.series_color(1),
            label="uniform plate of the same mass")
    for _, row in arm.iterrows():
        ax.annotate(
            f"{row['vf']:.2f}", (row["mass_kg"], row["compliance_J"]),
            textcoords="offset points", xytext=(0, 7), ha="center", fontsize=7.4,
            color=st.INK_SECONDARY,
        )
    ax.set_yscale("log")
    ax.set_ylabel("weighted compliance [J]")
    st.title(
        ax, "Aerospace bracket: mass-stiffness trade",
        "points are labelled with the volume-fraction target; lower is stiffer",
    )
    st.legend(ax, loc="upper right")

    ax = axes[1]
    if "f1_topology_Hz" in arm and arm["f1_topology_Hz"].notna().any():
        ax.plot(arm["mass_kg"], arm["f1_topology_Hz"], "o-",
                color=st.series_color(2),
                label="first natural frequency of the interpreted topology")
        if arm["f1_solid_Hz"].notna().any():
            ax.axhline(float(arm["f1_solid_Hz"].dropna().iloc[0]),
                       color=st.INK_MUTED, linewidth=1.0, linestyle="--",
                       label="full solid domain (also the equal-mass uniform plate)")
        ax.set_ylabel("f1 [Hz]")
    ax.set_xlabel("mass [kg]")
    st.title(
        ax, "vibration metric across the same sweep",
        "in this 2-D idealisation a uniformly thinned plate keeps the full-solid "
        "frequencies exactly, so the dashed line is the equal-mass baseline at "
        "every mass",
    )
    st.legend(ax, loc="lower right")
    st.annotate_note(
        fig,
        "The equal-mass uniform plate is C_solid / volume_fraction, which follows "
        "from compliance scaling as 1/thickness. Both curves come from solver "
        "runs recorded in the study summaries; no fit or extrapolation is used.",
    )
    return st.save_figure(fig, path)


def plot_numerical_settings(table: pd.DataFrame, path: str) -> str:
    """How the numerical settings change compliance and the design."""
    arms = [
        ("penalty", "SIMP penalty p [-]", "penalty"),
        ("filter_radius", "filter radius [cells]", "tag"),
        ("mesh_fixed_r", "elements", "num_elements"),
        ("youngs_modulus", "Young's modulus [Pa]", "youngs_modulus_Pa"),
    ]
    present = [a for a in arms if not table[table["arm"] == a[0]].empty]
    fig, axes = st.figure(7.4, 2.3 * len(present), nrows=len(present), ncols=1)
    axes = np.atleast_1d(axes).ravel()

    for ax, (arm, xlabel, key) in zip(axes, present):
        sub = table[table["arm"] == arm].copy()
        sub["x"] = sub[key] if key != "tag" else sub["tag"].map(_numeric_tag)
        sub = sub.sort_values("x")
        ax.plot(sub["x"], sub["compliance_J"], "o-", color=st.series_color(0),
                label="weighted compliance [J]")
        ax.set_ylabel("compliance [J]")
        ax.set_xlabel(xlabel)
        if arm in ("mesh_fixed_r", "youngs_modulus"):
            ax.set_xscale("log")
        twin_label = "grey level [-]"
        # No second y-axis: the grey level is shown as a direct label per point
        # so the two quantities never share a scale.
        for _, row in sub.iterrows():
            ax.annotate(
                f"grey {row['grey_level']:.2f}", (row["x"], row["compliance_J"]),
                textcoords="offset points", xytext=(0, 8), ha="center",
                fontsize=7.0, color=st.INK_SECONDARY,
            )
        st.title(ax, f"effect of {xlabel.split('[')[0].strip()}",
                 f"point labels give the {twin_label} of the resulting design")
        st.legend(ax, loc="best")

    st.annotate_note(
        fig,
        "Each arm changes one setting against the shared baseline deck. "
        "Compliance across the Young's-modulus arm scales as 1/E exactly, which "
        "is the expected linear-elastic behaviour and a useful consistency check "
        "rather than a design result.",
    )
    return st.save_figure(fig, path)


def plot_mesh_dependence(table: pd.DataFrame, path: str) -> str:
    """Mesh refinement with the filter radius fixed in metres vs in cells."""
    fixed_r = table[table["arm"] == "mesh_fixed_r"].copy()
    fixed_cells = table[table["arm"] == "mesh_fixed_cells"].copy()
    if fixed_r.empty or fixed_cells.empty:
        raise ValueError("the study table is missing a mesh-refinement arm")

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.6)

    ax = axes[0]
    for slot, (sub, label) in enumerate(
        [
            (fixed_r, "filter radius fixed at 3.75 mm (length scale held)"),
            (fixed_cells, "filter radius fixed at 1.5 cells (length scale shrinks)"),
        ]
    ):
        sub = sub.sort_values("num_elements")
        ax.plot(sub["num_elements"], sub["compliance_J"], "o-",
                color=st.series_color(slot), label=label)
    ax.set_xscale("log")
    ax.set_ylabel("compliance [J]")
    st.title(
        ax, "Mesh dependence of the optimised design",
        "holding the filter radius in physical units makes the optimum mesh "
        "convergent; holding it in cells does not",
    )
    st.legend(ax, loc="upper right")

    ax = axes[1]
    for slot, (sub, label) in enumerate(
        [(fixed_r, "radius fixed in metres"), (fixed_cells, "radius fixed in cells")]
    ):
        sub = sub.sort_values("num_elements")
        ax.plot(sub["num_elements"], sub["grey_level"], "o-",
                color=st.series_color(slot), label=label)
    ax.set_xscale("log")
    ax.set_xlabel("number of elements")
    ax.set_ylabel("grey level [-]")
    st.title(
        ax, "how binary the design is under refinement",
        "a radius fixed in cells produces ever-finer members, which is the "
        "classical mesh-dependence pathology filtering is meant to remove",
    )
    st.legend(ax, loc="best")
    return st.save_figure(fig, path)


def plot_load_weighting(table: pd.DataFrame, path: str) -> str:
    """Effect of the lateral load-case weight, objective and per-case.

    The sweep holds the down and reversal weights at (1, 0.5) and varies the
    lateral weight, so the x axis is the recorded weight of the lateral case
    rather than anything parsed out of a directory name. The one point that
    drops the reversal case entirely is a different experiment and is excluded
    from the curves; it is reported in the table and the topology montage.
    """
    arm = table[table["arm"] == "load_weighting"].copy()
    if arm.empty:
        raise ValueError("the study table has no load_weighting arm")

    required = ["weight_lateral", "weight_up_reversal", "compliance_down_limit_J",
                "compliance_up_reversal_J", "compliance_lateral_J"]
    missing = [c for c in required if c not in arm.columns]
    if missing:
        raise KeyError(
            "the load_weighting summaries carry no per-load-case data "
            f"(missing {missing}); re-run the study with a build that records "
            "optimization_setup.load_case_names"
        )

    sweep = arm[np.isclose(arm["weight_up_reversal"], 0.5)].sort_values("weight_lateral")
    if sweep.empty:
        raise ValueError("no load_weighting point holds the reversal weight at 0.5")

    fig, axes = st.stacked_panels(2, width=7.2, panel_height=2.6)

    ax = axes[0]
    ax.plot(sweep["weight_lateral"], sweep["compliance_J"], "o-",
            color=st.series_color(0),
            label="objective: weight-normalised mean over the three cases")
    last = len(sweep) - 1
    for position, (_, row) in enumerate(sweep.iterrows()):
        # End labels are aligned inwards so they cannot fall off the axes.
        align = "left" if position == 0 else "right" if position == last else "center"
        ax.annotate(
            f"grey {row['grey_level']:.3f}",
            (row["weight_lateral"], row["compliance_J"]),
            textcoords="offset points", xytext=(0, -13), ha=align, fontsize=7.0,
            color=st.INK_SECONDARY,
        )
    ax.set_ylabel("objective [J]")
    st.title(
        ax, "Aerospace bracket: effect of the lateral load-case weight",
        "weights are (down 1.0, reversal 0.5, lateral w). The objective is a "
        "mean over a weight set that changes along the axis, so it is not a "
        "like-for-like stiffness measure",
    )
    st.legend(ax, loc="upper right")

    # Each case relative to its value on the w = 0 design. The three absolute
    # levels differ by a factor of 19, so on a shared axis - log or not - the
    # changes that matter are invisible; the ratio is what the arm is about,
    # and the absolute values are in the table beside this figure.
    ax = axes[1]
    st.require_scatter_series(3, "load cases")
    reference = sweep[np.isclose(sweep["weight_lateral"], 0.0)]
    # The reversal case is exactly -0.5 times the down case, so their ratio
    # curves coincide; a dashed overlay keeps both readable instead of one
    # hiding the other.
    for slot, (column, label, style) in enumerate(
        [
            ("compliance_down_limit_J", "down limit case (weight 1.0)", "o-"),
            ("compliance_up_reversal_J",
             "up reversal case (weight 0.5) - coincides exactly", "o--"),
            ("compliance_lateral_J", "lateral case (weight w)", "o-"),
        ]
    ):
        base = float(reference[column].iloc[0]) if not reference.empty else 1.0
        ax.plot(sweep["weight_lateral"], sweep[column] / base, style,
                color=st.series_color(slot), label=label,
                markersize=7.5 if slot == 0 else 4.5,
                markerfacecolor="none" if slot == 1 else None,
                markeredgewidth=1.4 if slot == 1 else 0.0)
    ax.axhline(1.0, color=st.INK_MUTED, linewidth=0.9, linestyle="--")
    ax.set_xlabel("lateral load-case weight w [-]")
    ax.set_ylabel("compliance / its value at w = 0 [-]")
    st.title(
        ax, "what the weighting actually buys",
        "each series is one load case solved on the design that came out of "
        "that run. The lateral case gains 39 %; the two vertical cases do not "
        "pay for it, which is the local-optimum result discussed in the text",
    )
    st.legend(ax, loc="lower left")

    other = arm[~np.isclose(arm["weight_up_reversal"], 0.5)]
    note = (
        "Weights are normalised to sum to one inside the objective, so the "
        "objective is a mean and is not comparable in magnitude across "
        "different weight sets. The lower panel is each case relative to its "
        "own value on the w = 0 design; absolute compliances are tabulated in "
        "docs/aerospace_study.md."
    )
    if not other.empty:
        tags = ", ".join(
            f"{row['tag']} (objective {row['compliance_J']:.4g} J)"
            for _, row in other.iterrows()
        )
        note += (
            " Excluded from the curves because it changes a second weight: "
            f"{tags}."
        )
    st.annotate_note(fig, note)
    return st.save_figure(fig, path)


def plot_study_topologies(study_dir: str, table: pd.DataFrame, arm: str, path: str,
                          title_text: str, subtitle: str,
                          order_by: Optional[str] = None) -> str:
    """Small multiples of the optimised density for one study arm.

    `order_by` names the recorded column the panels are ordered by - the
    quantity the arm actually varies. Reading the order out of the directory
    names is the fallback.
    """
    from . import fields as fld
    from .loaders import load_mesh

    sub = table[table["arm"] == arm].copy()
    if sub.empty:
        raise ValueError(f"the study table has no {arm} arm")
    if order_by is not None and order_by in sub.columns and sub[order_by].notna().all():
        sub["x"] = sub[order_by]
    else:
        sub["x"] = sub["tag"].map(_numeric_tag)
    sub = sub.sort_values(["x", "tag"])

    count = len(sub)
    columns = min(4, count)
    rows = (count + columns - 1) // columns
    fig, axes = st.figure(7.6, 1.6 * rows + 0.4, nrows=rows, ncols=columns)
    axes = np.atleast_1d(axes).ravel()

    for slot, (_, row) in enumerate(sub.iterrows()):
        directory = os.path.join(study_dir, arm, row["tag"])
        mesh = load_mesh(directory)
        density = load_csv(os.path.join(directory, "density_final.csv"))
        values = density["physical_density[-]"].to_numpy()
        ax = axes[slot]
        shape = fld.structured_grid_shape(mesh)
        extent = mesh.extent
        if shape is not None:
            ax.imshow(fld.density_image(values, shape), cmap=st.DENSITY_CMAP_SURFACE,
                      vmin=0.0, vmax=1.0, origin="lower",
                      extent=(extent[0], extent[1], extent[2], extent[3]),
                      interpolation="nearest", aspect="equal")
        else:
            fld.element_collection(ax, mesh, values, cmap=st.DENSITY_CMAP_SURFACE,
                                   vmin=0.0, vmax=1.0)
            fld.set_domain_limits(ax, mesh)
        fld.bare_axes(ax)
        ax.set_title(
            f"{row['tag']}\nc = {st.format_si(row['compliance_J'])} J",
            loc="left", fontsize=8.4,
        )
    for slot in range(count, len(axes)):
        axes[slot].axis("off")

    fig.suptitle(title_text, x=0.01, ha="left", fontsize=11, fontweight="bold",
                 color=st.INK_PRIMARY)
    st.annotate_note(fig, subtitle)
    return st.save_figure(fig, path)


# ---------------------------------------------------------------------------
# Quadratic tetrahedra and linear buckling
# ---------------------------------------------------------------------------
#: Fixed categorical slot per element for the Tet10 and buckling figures, so
#: an element keeps its colour from one figure to the next.
ELEMENT_SLOTS = {"Hex8": 0, "Tet4": 1, "Tet10": 2, "Quad4": 3}


def plot_tet10_convergence(directory: str, path: str) -> str:
    """Cantilever tip error of Hex8, Tet4 and Tet10 on the same grids."""
    table = load_csv(os.path.join(directory, "mesh_convergence_tet10.csv"))
    summary = load_json(os.path.join(directory, "summary.json")).get(
        "mesh_convergence_tet10", {})
    fig, ax = st.figure(7.4, 4.4)
    for element in ("Hex8", "Tet4", "Tet10"):
        sub = table[table["element"] == element].sort_values("num_dofs")
        if sub.empty:
            continue
        ax.plot(sub["num_dofs"], sub["rel_error_timoshenko[-]"],
                ELEMENT_MARKERS[element] + "-", color=st.series_color(ELEMENT_SLOTS[element]),
                label=element)
        last = sub.iloc[-1]
        ax.annotate(f"{100.0 * float(last['rel_error_timoshenko[-]']):.2g} %",
                    (float(last["num_dofs"]), float(last["rel_error_timoshenko[-]"])),
                    xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                    color=st.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("relative error vs Timoshenko [-]")
    st.title(ax, "Tip deflection error: linear against quadratic elements",
             "solid cantilever 1 x 0.1 x 0.05 m, nu = 0, clamped root, tip traction with "
             "a 1 kN resultant; grids 10 x 2 x 1, 20 x 4 x 2, 40 x 8 x 4, each cell one "
             "Hex8 or six Kuhn tetrahedra")
    st.legend(ax, loc="lower left")
    st.annotate_note(
        fig,
        "The quadratic tetrahedron contains the linear-stress field of pure bending, so "
        "it does not lock: on the coarsest grid it is already within "
        f"{100.0 * float(summary.get('Tet10', {}).get('records', [{}])[0].get('relative_error_vs_timoshenko', float('nan'))):.2g} % "
        "of beam theory. The error is measured against Timoshenko beam theory, so it "
        "also holds the modelling gap between a beam and solid elasticity.",
    )
    return st.save_figure(fig, path)


def plot_buckling_verification(directory: str, path: str) -> str:
    """Critical load of a clamped column against Euler-Engesser."""
    table = load_csv(os.path.join(directory, "buckling_euler.csv"))
    fig, ax = st.figure(7.4, 4.4)
    for element in ("Quad4", "Hex8", "Tet4", "Tet10"):
        sub = table[table["element"] == element].sort_values("num_dofs")
        if sub.empty:
            continue
        label = "Q4 (plane stress)" if element == "Quad4" else element
        ax.plot(sub["num_dofs"], sub["rel_error_engesser[-]"],
                ELEMENT_MARKERS[element] + "-", color=st.series_color(ELEMENT_SLOTS[element]),
                label=label)
        last = sub.iloc[-1]
        ax.annotate(f"{100.0 * float(last['rel_error_engesser[-]']):.2g} %",
                    (float(last["num_dofs"]), float(last["rel_error_engesser[-]"])),
                    xytext=(6, 0), textcoords="offset points", va="center", fontsize=8,
                    color=st.INK_SECONDARY)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("(lambda_1 - P_Engesser) / P_Engesser [-]")
    # The Tet10 meshes halve h exactly, so Richardson extrapolation of the
    # three errors separates the discretisation error from the modelling gap.
    gap_text = ""
    tet10 = table[table["element"] == "Tet10"].sort_values("num_dofs")
    if len(tet10) >= 3:
        e1, e2, e3 = tet10["rel_error_engesser[-]"].to_numpy()[-3:]
        if (e1 - e2) * (e2 - e3) > 0 and e2 != e3:
            ratio = (e1 - e2) / (e2 - e3)
            if ratio > 1.0:
                order = np.log2(ratio)
                limit = e3 - (e2 - e3) / (ratio - 1.0)
                gap_text = (f" Richardson extrapolation of the Tet10 errors (order "
                            f"{order:.2f} in h) puts their limit at {100.0 * limit:.2f} %: the "
                            "gap between the solid model and the beam formula, not a "
                            "discretisation error.")
    st.title(ax, "Linear buckling of a clamped column: first load factor vs mesh",
             "1 m steel column, 50 mm square section (3-D) or 50 x 20 mm (plane stress), "
             "unit axial tip traction, so lambda_1 is the critical load in newtons; "
             "reference Euler P = pi^2 E I / (4 L^2) with Engesser's shear correction")
    st.legend(ax, loc="lower left")
    st.annotate_note(
        fig,
        "Every point lies above the reference: a displacement-based element is too "
        "stiff, and its buckling load converges from above. The linear Hex8 and Tet4 "
        "lock in bending and converge slowly; the Q4 and the quadratic Tet10 reach "
        "the 1 % tolerance of sparlab_verify's buckling-euler study." + gap_text,
    )
    return st.save_figure(fig, path)


def plot_tet10_part_study(directory: str, path: str) -> str:
    """Compliance of the engine mount against DOFs, per element and geometry."""
    table = load_csv(os.path.join(directory, "tet10_part_study.csv"))
    summary = load_json(os.path.join(directory, "summary.json"))
    variants = [("Tet4", "^", 1), ("Tet10 (straight)", "d", 2), ("Tet10 (curved)", "D", 3)]
    fig, axes = st.figure(7.4, 6.6, nrows=2, ncols=1, sharex=True)
    for ax, case, label in zip(axes, ("vertical", "lateral"),
                               ("vertical pin load, 12 kN", "lateral pin load, 5 kN")):
        key = f"compliance_{case}_J"
        for variant, marker, slot in variants:
            sub = table[table["variant"] == variant].sort_values("num_dofs")
            if sub.empty:
                continue
            ax.plot(sub["num_dofs"], sub[key], marker + "-", color=st.series_color(slot - 1),
                    label=variant)
            for _, row in sub.iterrows():
                ax.annotate(f"{row['size_mm']:g}", (row["num_dofs"], row[key]),
                            xytext=(0, 6 if variant == "Tet10 (curved)" else -11),
                            textcoords="offset points", ha="center", fontsize=7,
                            color=st.INK_MUTED)
        ref = summary["references"][case]
        ax.axhline(ref["finest_tet10_curved_J"], color=st.INK_MUTED, linewidth=0.9,
                   linestyle="--")
        ax.set_xscale("log")
        ax.set_ylabel("compliance f.u [J]")
        st.title(ax, f"Engine mount, {label}",
                 f"dashed: the finest curved Tet10 run ({ref['finest_size_mm']:g} mm), itself "
                 "a lower bound; numbers at the points are the element size in mm")
        st.legend(ax, loc="lower right")
    axes[-1].set_xlabel("number of degrees of freedom")
    st.annotate_note(
        fig,
        "Each mesh comes from the same CAD model through Gmsh with the curvature "
        "refinement scaled with the element size. A displacement-based solution is "
        "too stiff, so its compliance lies below the exact value and rises towards it. "
        "The straight-sided Tet10 meshes are the linear meshes with an edge node at "
        "every edge midpoint: the same faceted holes as the Tet4 mesh, a slightly "
        "different part from the curved one.",
    )
    return st.save_figure(fig, path)


def plot_buckling_cross_validation(directory: str, path: str) -> str:
    """Buckling load factors against scikit-fem and CalculiX *BUCKLE."""
    summary = load_json(os.path.join(directory, "summary.json"))
    codes = [("scikit-fem buckling", "scikit-fem (same K_G, dense eigensolve)"),
             ("calculix buckling", "CalculiX *BUCKLE")]
    rows = []
    for case in summary.get("cases", []):
        for load_case in case.get("load_cases", []):
            if any(code in load_case["codes"] for code, _ in codes):
                rows.append((case["case"], case.get("element_type", ""),
                             load_case["load_case"], load_case["codes"]))
    if not rows:
        raise ValueError("cross-validation summary has no buckling comparisons")
    fig, ax = st.figure(7.4, 0.5 * len(rows) + 3.3)
    y = np.arange(len(rows))[::-1]
    for slot, (code, label) in enumerate(codes):
        xs, ys = [], []
        for position, (*_rest, results) in zip(y, rows):
            entry = results.get(code)
            if entry is not None:
                xs.append(entry["max_rel_diff"])
                ys.append(position)
        if xs:
            ax.plot(xs, ys, "o", color=st.series_color(slot), markersize=7, label=label)
    tolerances = summary.get("tolerances", {})
    for slot, key, label in ((0, "skfem_buckling", "scikit-fem tolerance"),
                             (1, "calculix_buckling", "CalculiX tolerance")):
        value = tolerances.get(key)
        if value:
            ax.axvline(value, color=st.series_color(slot), linewidth=1.0, linestyle="--",
                       label=f"{label} {value:g}")
    ax.set_xscale("log")
    ax.set_xlim(1.0e-12, 1.0e-2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c} ({e}), '{l}'" for c, e, l, _r in rows], fontsize=8.0)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("max over the reported modes of |lambda_ref - lambda| / lambda [-]")
    st.title(ax, "Cross-validation: linear buckling load factors",
             "scikit-fem assembles the geometric stiffness from its own static solution at "
             "the same quadrature points (the same discrete problem); CalculiX runs "
             "*BUCKLE on the exported deck with its own stress-stiffness evaluation")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=8,
              frameon=False)
    st.annotate_note(
        fig,
        "Where an element's stress varies inside it (C3D8, C3D10), CalculiX's factors "
        "lie a few 1e-5 above SparLab's and scikit-fem's, which agree to about 1e-9; "
        "for C3D4, whose stress is constant in an element, the gap is several times "
        "smaller. Which detail of CalculiX's stress stiffness accounts for it has not "
        "been identified; the tolerance records the measured size, it does not "
        "explain it.",
    )
    return st.save_figure(fig, path)
