"""Figures for topology optimisation under loads that follow the design:
self-weight, and heating of a structure held at its ends.

Every number on these figures is read from a run's `summary.json` or its
density file; nothing is recomputed. Each figure sets one deck beside the
comparison runs `scripts/run_all_benchmarks.sh` makes of it with the load
scaled (`--gravity-scale`, `--temperature-scale`).
"""

from __future__ import annotations

from typing import Optional, Sequence, Tuple

import numpy as np

from . import fields as fld
from . import style as st
from .buckling import part_load_factor
from .loaders import CaseResults


def _density(case: CaseResults) -> np.ndarray:
    table = case.density()
    if table is None:
        raise FileNotFoundError(f"{case.directory} has no density_final.csv")
    return table["physical_density[-]"].to_numpy()


def _result(case: CaseResults) -> dict:
    return case.summary.get("optimization_result", {})


def _target(case: CaseResults) -> float:
    return float(case.summary.get("optimization_setup", {}).get("volume_fraction_target",
                                                                float("nan")))


def _design_panel(ax, case: CaseResults, title: str) -> None:
    density = _density(case)
    fld.element_collection(ax, case.mesh, density, cmap=st.DENSITY_CMAP_SURFACE,
                           vmin=0.0, vmax=1.0)
    fld.mesh_outline(ax, case.mesh, color=st.INK_MUTED, linewidth=0.6)
    fld.set_domain_limits(ax, case.mesh, margin=0.03)
    fld.bare_axes(ax)
    ax.set_title(title, loc="left", fontsize=8.6)


def plot_self_weight_comparison(designs: Sequence[Tuple[float, CaseResults]], path: str,
                                name: Optional[str] = None) -> str:
    """One deck at several multiples of gravity: the designs side by side, and
    the volume fraction each uses of its allowance."""
    if not designs:
        raise ValueError("no designs to compare")
    for _, case in designs:
        if case.dim != 2:
            raise ValueError("the self-weight figure is drawn for plane cases")
    count = len(designs)
    xmin, xmax, ymin, ymax = designs[0][1].mesh.extent
    aspect = (ymax - ymin) / max(xmax - xmin, 1.0e-12)
    field = float(np.clip(7.4 / count * aspect * 0.9, 1.0, 3.2))
    fig, axes = st.figure(7.4, field + 4.2, nrows=2, ncols=count,
                          height_ratios=[field, 2.4])
    axes = np.asarray(axes).reshape(2, count)
    for ax in axes[1, 1:]:
        ax.remove()
    axes[1, 0].remove()
    chart = fig.add_subplot(axes[0, 0].get_gridspec()[1, :])

    used, targets, labels = [], [], []
    for k, (scale, case) in enumerate(designs):
        result = _result(case)
        label = "no self-weight" if scale == 0 else f"{scale:g} g"
        _design_panel(axes[0, k], case,
                      f"{label}\ncompliance {st.format_si(result.get('compliance_J', 0.0))} J")
        used.append(float(result.get("volume_fraction", float("nan"))))
        targets.append(_target(case))
        labels.append(label)

    x = np.arange(count)
    chart.bar(x, used, width=0.5, color=st.series_color(0), label="volume fraction used")
    for i, v in enumerate(used):
        chart.annotate(f"{v:.3f}", (x[i], v), xytext=(0, 3), textcoords="offset points",
                       ha="center", va="bottom", fontsize=8, color=st.INK_SECONDARY)
    target = targets[0]
    chart.axhline(target, color=st.INK_SECONDARY, linewidth=1.0, linestyle="--")
    chart.annotate(f"allowed {target:g}", (count - 0.5, target), xytext=(0, 4),
                   textcoords="offset points", ha="right", va="bottom", fontsize=8,
                   color=st.INK_SECONDARY)
    chart.set_xticks(x)
    chart.set_xticklabels(labels, fontsize=8)
    chart.set_xlim(-0.6, count - 0.4)
    chart.set_ylim(0.0, 1.3 * max(target, max(used)))
    chart.set_ylabel("volume fraction [-]")
    st.title(chart, "Material used",
             "the weight grows with the material it moves, so under enough of it the "
             "optimum stops using its allowance")

    first = designs[0][1]
    st.figure_title(
        fig, f"{name or first.name}: minimum compliance under a deck load and a weight "
        "that follows the design",
        "same mesh, supports, deck load, filter, projection and allowance; gravity "
        "scaled from 0 (the deck load alone) upwards; physical density after the "
        "projection",
    )
    st.annotate_note(
        fig,
        "The weight of each cell is gamma(rho) times the solid's, gamma = rho above the "
        "body-load threshold (0.1) and falling like rho^p below it. The compliances of "
        "different loads are not comparable with one another; each is the objective "
        "its own run minimised.",
    )
    return st.save_figure(fig, path)


def plot_thermal_sweep(designs: Sequence[Tuple[float, CaseResults]], path: str,
                       shown: Sequence[float] = (0.0, 1.0, 3.0, 5.0, 10.0, 40.0),
                       name: Optional[str] = None) -> str:
    """One deck at a range of temperature rises: the designs at `shown`, the
    volume fraction each uses and the lowest buckling load factor of each
    exported part."""
    if not designs:
        raise ValueError("no designs to compare")
    designs = sorted(designs, key=lambda d: d[0])
    for _, case in designs:
        if case.dim != 2:
            raise ValueError("the thermal figure is drawn for plane cases")
    picked = [d for d in designs if any(abs(d[0] - s) < 1.0e-9 for s in shown)]
    if not picked:
        raise ValueError("none of the shown temperature rises was run")
    columns = 3
    rows = (len(picked) + columns - 1) // columns
    xmin, xmax, ymin, ymax = designs[0][1].mesh.extent
    aspect = (ymax - ymin) / max(xmax - xmin, 1.0e-12)
    field = float(np.clip(7.4 / columns * aspect * 1.1, 0.8, 2.6))
    fig, axes = st.figure(7.4, rows * (field + 0.45) + 3.9, nrows=rows + 1, ncols=columns,
                          height_ratios=[field] * rows + [2.5])
    axes = np.asarray(axes).reshape(rows + 1, columns)
    for k in range(rows * columns):
        ax = axes[k // columns, k % columns]
        if k >= len(picked):
            ax.remove()
            continue
        rise, case = picked[k]
        result = _result(case)
        _design_panel(ax, case,
                      f"dT = {rise:g} K\nvolume fraction "
                      f"{float(result.get('volume_fraction', float('nan'))):.3f}")
    for ax in axes[rows, :]:
        ax.remove()
    grid = axes[0, 0].get_gridspec()
    left = fig.add_subplot(grid[rows, 0:2])
    right = fig.add_subplot(grid[rows, 2])

    rises = np.array([d[0] for d in designs])
    used = np.array([float(_result(c).get("volume_fraction", float("nan")))
                     for _, c in designs])
    lam = np.array([part_load_factor(c) or float("nan") for _, c in designs])
    target = _target(designs[0][1])
    left.plot(rises, used, "o-", color=st.series_color(0), markersize=5,
              label="volume fraction used")
    left.axhline(target, color=st.INK_SECONDARY, linewidth=1.0, linestyle="--")
    left.annotate(f"allowed {target:g}", (rises.max(), target), xytext=(0, 4),
                  textcoords="offset points", ha="right", va="bottom", fontsize=8,
                  color=st.INK_SECONDARY)
    left.set_xscale("symlog", linthresh=1.0)
    left.set_xlim(-0.2, rises.max() * 1.15)
    left.set_ylim(0.0, 1.25 * target)
    left.set_xlabel("temperature rise dT [K]")
    left.set_ylabel("volume fraction [-]")
    # Where the volume constraint stops binding, as the runs show it.
    active = used >= target * (1.0 - 1.0e-2)
    turn = next((i for i in range(1, len(rises)) if active[i - 1] and not active[i]), None)
    st.title(left, "Material used",
             (f"the volume constraint turns inactive between {rises[turn - 1]:g} and "
              f"{rises[turn]:g} K" if turn is not None else
              "the volume fraction each design uses of its allowance"),
             wrap=60)
    right.plot(rises, lam, "s-", color=st.series_color(1), markersize=5)
    right.set_xscale("symlog", linthresh=1.0)
    right.set_yscale("log")
    right.set_xlim(-0.2, rises.max() * 1.15)
    right.set_xlabel("dT [K]")
    right.set_ylabel("lowest load factor [-]")
    st.title(right, "Exported part", "linear buckling", wrap=30)

    first = designs[0][1]
    st.figure_title(
        fig, f"{name or first.name}: minimum compliance of a beam held at both ends, "
        "loaded at its centre and heated",
        "same mesh, supports, load, filter, projection and allowance; the temperature "
        "rise swept; physical density after the projection",
    )
    st.annotate_note(
        fig,
        "The thermal load of each cell is E(rho)/E0 times the solid's, and the objective "
        "is the compliance of the whole load, mechanical and thermal. A statically "
        "determinate truss lets its members expand freely and so escapes the thermal "
        "part; the thinner members it keeps buckle at a smaller multiple of the design "
        "load.",
    )
    return st.save_figure(fig, path)
