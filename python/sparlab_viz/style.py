"""Figure style: one place for colour, type and layout decisions.

Colour choices follow the job the colour does, not taste:

* **Categorical** (which series) - the eight-slot order below, assigned in
  sequence and never cycled. It is validated for lightness band, chroma floor,
  protan/deuteran separation and normal-vision separation on a light surface
  (worst adjacent pair Delta-E 9.1 CVD / 19.6 normal in OKLab x100). Three
  slots also clear the harder all-pairs test, which is why scatter-style
  figures here stay at three series or fewer.
* **Sequential** (how much: displacement magnitude, von Mises, strain energy) -
  `viridis`. Monotone in lightness, perceptually uniform and readable under
  colour-vision deficiency; never a rainbow map.
* **Density** - `Greys`, so 1.0 prints black. Achromatic is the strongest
  possible single-hue ramp and matches the field convention for SIMP designs.
* **Diverging** (signed stress components, sensitivity error sign) - `coolwarm`,
  which has two hues either side of a neutral grey midpoint. Limits are always
  made symmetric about zero so the midpoint really is zero.

Two further rules are enforced by construction rather than by review:

* no chart in this package uses two y-axes. Quantities of different scale go in
  separate stacked panels sharing the x axis;
* every figure with two or more series carries a legend, so identity is never
  colour-alone.
"""

from __future__ import annotations

import os
import textwrap
from typing import Optional, Sequence

import matplotlib
matplotlib.use("Agg")  # figures are written to files, never shown

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize, TwoSlopeNorm

# --- surfaces and ink -------------------------------------------------------
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8982"
GRID = "#e3e2dc"

# --- categorical slots, fixed order ----------------------------------------
SERIES = [
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#4a3aa7",  # 7 violet
    "#e34948",  # 8 red
]

#: Number of series a scatter-style figure may carry: only the first three
#: categorical slots clear the all-pairs separation test.
MAX_SCATTER_SERIES = 3

# --- ramps ------------------------------------------------------------------
FIELD_CMAP = "viridis"
SIGNED_CMAP = "coolwarm"
DENSITY_CMAP = "Greys"

#: A density ramp that keeps void areas on the page surface rather than pure
#: white, so a figure reads correctly when printed on off-white paper.
DENSITY_CMAP_SURFACE = LinearSegmentedColormap.from_list(
    "sparlab_density", [SURFACE, "#111111"]
)


def series_color(index: int) -> str:
    """Categorical colour for series `index`, assigned in fixed order.

    A ninth series is a modelling decision, not a colour decision: fold it into
    an "other" group or split the figure into small multiples.
    """
    if index >= len(SERIES):
        raise ValueError(
            f"series index {index} exceeds the {len(SERIES)}-slot categorical "
            "palette; fold extra series together or use small multiples instead "
            "of generating a new hue"
        )
    return SERIES[index]


def apply_style() -> None:
    """Install the SparLab rcParams. Idempotent."""
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "figure.dpi": 130,
            "savefig.dpi": 160,
            "savefig.facecolor": SURFACE,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "axes.facecolor": SURFACE,
            "axes.edgecolor": INK_MUTED,
            "axes.linewidth": 0.8,
            "axes.labelcolor": INK_SECONDARY,
            "axes.titlecolor": INK_PRIMARY,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.titlepad": 8,
            "axes.labelsize": 9.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "axes.axisbelow": True,
            "grid.color": GRID,
            "grid.linewidth": 0.7,
            "xtick.color": INK_SECONDARY,
            "ytick.color": INK_SECONDARY,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "legend.frameon": False,
            "legend.fontsize": 8.5,
            "legend.labelcolor": INK_SECONDARY,
            "legend.handlelength": 1.8,
            "lines.linewidth": 1.9,
            "lines.markersize": 5.0,
            "lines.markeredgewidth": 0.0,
            "font.size": 9.5,
            "font.family": "sans-serif",
            "text.color": INK_PRIMARY,
            "image.cmap": FIELD_CMAP,
        }
    )


#: Characters per line used when wrapping subtitles and notes. Matplotlib's own
#: `wrap=True` measures against the *figure* width and ignores the axes, which
#: on a multi-panel figure lets a subtitle run under a colour bar; wrapping the
#: string ourselves is predictable.
SUBTITLE_WRAP = 96
NOTE_WRAP = 118


def figure(width: float = 7.0, height: float = 4.0, **kwargs):
    """Create a styled figure and axes.

    Constrained layout is used everywhere: these figures carry titles,
    subtitles, colour bars and footnotes, and hand-tuned spacing does not
    survive a change of panel count or domain aspect ratio.
    """
    apply_style()
    kwargs.setdefault("layout", "constrained")
    return plt.subplots(figsize=(width, height), **kwargs)


def stacked_panels(n: int, width: float = 7.0, panel_height: float = 2.2):
    """`n` panels sharing the x axis.

    This is the replacement for a two-scale (twinx) chart: quantities whose
    magnitudes differ get their own panel and their own y label, so no reader
    has to work out which curve belongs to which axis.
    """
    apply_style()
    fig, axes = plt.subplots(
        n, 1, figsize=(width, panel_height * n), sharex=True, layout="constrained"
    )
    if n == 1:
        axes = [axes]
    return fig, list(axes)


def title(ax, text: str, subtitle: Optional[str] = None) -> None:
    """Left-aligned panel title with an optional smaller subtitle beneath it.

    The subtitle is placed in offset *points* above the axes and the title is
    padded to clear it, so the two never overlap however tall the panel is.
    """
    if not subtitle:
        ax.set_title(text, loc="left")
        return
    wrapped = textwrap.fill(subtitle, SUBTITLE_WRAP)
    lines = wrapped.count("\n") + 1
    ax.annotate(
        wrapped, xy=(0.0, 1.0), xycoords="axes fraction", xytext=(0.0, 4.0),
        textcoords="offset points", fontsize=8.5, color=INK_SECONDARY,
        va="bottom", ha="left", annotation_clip=False,
    )
    ax.set_title(text, loc="left", pad=8.0 + 10.5 * lines)


def figure_title(fig, text: str, subtitle: Optional[str] = None, ax=None) -> None:
    """Title and subtitle as a reserved band across the top of the figure.

    Both lines are figure-level text and the layout engine is told to keep the
    panels below them (`rect`). Neither is an axes title or annotation, because
    on a multi-column grid a wide subtitle attached to one panel counts as that
    panel's decoration: the engine then shrinks every panel to fit the text
    inside one cell, which on an equal-aspect field panel collapses the figure.

    `ax` is accepted for call compatibility and is not used.
    """
    del ax  # the band is figure-level; no panel owns it
    height = float(fig.get_size_inches()[1])
    top_in = 0.07          # figure top to the top of the title line
    title_in = 0.21        # title line height at 11 pt bold
    line_in = 0.145        # subtitle line height at 8.5 pt
    gap_in = 0.10          # band to first panel

    band = top_in + title_in
    wrapped = textwrap.fill(subtitle, SUBTITLE_WRAP) if subtitle else ""
    if wrapped:
        band += line_in * (wrapped.count("\n") + 1)

    engine = fig.get_layout_engine()
    if engine is not None and hasattr(engine, "set"):
        engine.set(rect=(0.0, 0.0, 1.0, max(0.3, 1.0 - (band + gap_in) / height)))

    fig.text(0.012, 1.0 - top_in / height, text, ha="left", va="top",
             fontsize=11, fontweight="bold", color=INK_PRIMARY)
    if wrapped:
        fig.text(0.012 + 0.055, 1.0 - (top_in + title_in) / height, wrapped,
                 ha="left", va="top", fontsize=8.5, color=INK_SECONDARY)


def annotate_note(fig, text: str) -> None:
    """Add a footnote line, used for the interpretation caveats."""
    fig.supxlabel(
        textwrap.fill(text, NOTE_WRAP), fontsize=7.6, color=INK_MUTED,
        ha="left", x=0.004,
    )


def limit_ticks(ax, x: int = 5, y: int = 5) -> None:
    """Cap the tick count so labels cannot collide on a small panel."""
    from matplotlib.ticker import MaxNLocator

    ax.xaxis.set_major_locator(MaxNLocator(nbins=x, prune="both"))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=y))


def symmetric_norm(values, cmap_center: float = 0.0) -> Normalize:
    """Diverging normalisation with `cmap_center` exactly at the midpoint."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return Normalize(vmin=-1.0, vmax=1.0)
    span = max(abs(float(finite.min()) - cmap_center),
               abs(float(finite.max()) - cmap_center))
    if span == 0.0:
        span = 1.0
    return TwoSlopeNorm(vmin=cmap_center - span, vcenter=cmap_center,
                        vmax=cmap_center + span)


def save_figure(fig, path: str, close: bool = True) -> str:
    """Write a figure, creating its directory. Returns the path."""
    directory = os.path.dirname(os.path.abspath(path))
    if directory:
        os.makedirs(directory, exist_ok=True)
    fig.savefig(path)
    if close:
        plt.close(fig)
    return path


def format_si(value: float, digits: int = 4) -> str:
    """Compact number formatting for direct labels.

    Significant figures are kept including trailing zeros, because these
    labels sit beside the result tables: a compliance of 3.00015 J has to read
    `3.000`, not `3`, or the figure looks like it disagrees with the table.
    """
    if value == 0:
        return "0"
    magnitude = abs(value)
    if 1e-3 <= magnitude < 1e5:
        text = f"{value:#.{digits}g}"
        return text[:-1] if text.endswith(".") else text
    return f"{value:.{digits - 1}e}"


def require_scatter_series(count: int, what: str = "series") -> None:
    """Guard the three-series cap that all-pairs colour separation imposes."""
    if count > MAX_SCATTER_SERIES:
        raise ValueError(
            f"{count} {what} in a scatter-style figure exceeds the "
            f"{MAX_SCATTER_SERIES}-series cap that all-pairs colour separation "
            "allows; use small multiples or fold series together"
        )


def legend(ax, **kwargs) -> None:
    """Attach a legend with the house settings."""
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) >= 2:
        ax.legend(handles, labels, **kwargs)


def sequence_colors(n: int, cmap: str = "viridis", lo: float = 0.12,
                    hi: float = 0.92) -> Sequence:
    """Ordered colours for an *ordinal* series set (mesh sizes, iterations).

    Ordinal data takes a single-hue-style ramp so the reader sees the order in
    the colour; using categorical slots there would spend the identity channel
    on something the ordering already encodes.
    """
    return plt.get_cmap(cmap)(np.linspace(lo, hi, max(n, 1)))
