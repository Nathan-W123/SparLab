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
def plot_runtime_scaling(directory: str, path: str) -> str:
    """Wall-clock cost of each phase against problem size."""
    table = load_csv(os.path.join(directory, "runtime_scaling.csv"))
    summary = load_json(os.path.join(directory, "runtime_scaling.json"))
    exponents = summary.get("scaling_exponent_vs_dofs", {})

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
        ax, "Runtime scaling of the solver phases",
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
        "flat with size, as expected for a fixed-stencil structured Q4 mesh",
    )
    st.legend(ax, loc="lower right")

    st.annotate_note(
        fig,
        "Slopes are least-squares fits of log(time) against log(DOFs) over the "
        "largest three sizes. Assembly is O(n); a 2-D sparse Cholesky with a "
        "fill-reducing ordering is close to O(n^1.5) asymptotically, and at "
        "these sizes the measured slope also carries cache effects.",
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
