#!/usr/bin/env python3
"""Social-post cover image of an optimised part, from a real run.

usage:
    python3 python/scripts/make_cover_image.py \
        --case results/aerospace_bracket --load-case down_limit \
        --hook "What material does this structure actually need?" \
        --output docs/figures/cover_aerospace_bracket.png

Same data as `make_hero_figure.py` and none of the apparatus. The part fills
the frame on a dark surface, carrying its real stress field, with only the
supports, the applied load and one line of text. Axes, grid, legend, the design
domain and the explanatory note are all dropped: on a phone the structure has
to be the whole image.

Honesty notes, since a cover has no room for a caption:

* the field is `solid_von_mises`, the von Mises stress of the solid material,
  the same column the technical figure uses;
* the colour range is clipped to a percentile window of the retained material
  (2nd to 99th by default) rather than to 0 and the peak. On a dark surface a
  range starting at zero puts the low-stress material at viridis's near-black
  end, where it disappears into the background, and the true peak is a
  mesh-sensitive re-entrant-corner value that would flatten everything else.
  Both limits are printed on the colour key and reported on stdout;
* no claim is made by the picture. `--hook` text is the author's to get right;
  the one thing the run does support directly is how much material was removed,
  which `--hook-auto` writes for you.
"""

from __future__ import annotations

import argparse
import sys
import textwrap

import _bootstrap  # noqa: F401

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from make_hero_figure import _clusters, load_resultant, read_case

from sparlab_viz import fields as fld
from sparlab_viz import style as st
from sparlab_viz.loaders import ResultError

#: Dark surface. Not pure black: a near-black with a slight cool cast keeps the
#: low end of viridis distinguishable from the background.
COVER_SURFACE = "#0b0e13"
COVER_INK = "#f4f4f2"
COVER_INK_MUTED = "#8d93a1"



def _fit_text(fig, artist, max_fraction: float) -> float:
    """Shrink `artist` until it is no wider than `max_fraction` of the figure.

    Returns the height it ends up occupying, as a figure fraction, so the next
    line can be placed under it.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width_px = fig.get_size_inches()[0] * fig.dpi
    for _ in range(40):
        box = artist.get_window_extent(renderer)
        if box.width <= max_fraction * width_px:
            break
        artist.set_fontsize(artist.get_fontsize() * 0.94)
    box = artist.get_window_extent(renderer)
    return float(box.height / (fig.get_size_inches()[1] * fig.dpi))


def build(case_dir: str, load_case: str, path: str, hook: str = "",
          subhook: str = "", threshold: float = None,
          pixels: tuple = (1600, 1200), low_percentile: float = 2.0,
          high_percentile: float = 99.0, colorkey: bool = True,
          margin: float = 0.06) -> str:
    case, keep, values, threshold = read_case(case_dir, load_case, threshold)
    mesh = case.mesh
    if not keep.any():
        raise ResultError("no element reaches the interpretation threshold")

    shown = values[keep]
    vmin = float(np.percentile(shown, low_percentile))
    vmax = float(np.percentile(shown, high_percentile))

    # Frame the material, not the design domain: the empty part of the domain is
    # exactly what the cover should not spend area on.
    polygons = mesh.polygons[keep]
    xs, ys = polygons[:, :, 0], polygons[:, :, 1]
    xmin, xmax = float(xs.min()), float(xs.max())
    ymin, ymax = float(ys.min()), float(ys.max())

    width_px, height_px = pixels
    dpi = 200.0
    fig = plt.figure(figsize=(width_px / dpi, height_px / dpi), dpi=dpi,
                     facecolor=COVER_SURFACE)
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.set_facecolor(COVER_SURFACE)
    ax.set_axis_off()

    fld.element_collection(ax, mesh, values, cmap=st.FIELD_CMAP,
                           vmin=vmin, vmax=vmax, mask=keep)

    # --- supports and load, the two marks that carry the engineering story ---
    fixed = mesh.constrained_nodes()
    span = max(xmax - xmin, ymax - ymin)
    if fixed.size:
        for centre, radius in _clusters(mesh.nodes[fixed]):
            ax.add_patch(Circle(centre, radius, facecolor="none",
                                edgecolor=st.series_color(4), linewidth=3.4,
                                zorder=4))
    applied = load_resultant(mesh, load_case)
    if applied is not None:
        tail, force, _count = applied
        magnitude = float(np.hypot(*force))
        head = tail + force / max(magnitude, 1.0e-30) * 0.26 * span
        ax.annotate("", xy=tuple(head), xytext=tuple(tail),
                    annotation_clip=False,
                    arrowprops=dict(arrowstyle="-|>", linewidth=4.0,
                                    color=st.series_color(1), shrinkA=0,
                                    shrinkB=0, mutation_scale=34))

    # --- framing ----------------------------------------------------------
    # The structure is fitted to the frame, minus a band at the top reserved for
    # the hook, so it occupies as much of the image as its own aspect allows
    # instead of sitting in a letterboxed strip. Equal aspect means the view box
    # must match the frame's aspect, so the fit is done on the view height.
    frame_aspect = width_px / height_px
    band_top = 0.17 if hook else 0.0
    if subhook:
        band_top += 0.06
    # The colour key gets a band of its own rather than being dropped on top of
    # the structure, which fills everything the hook does not.
    band_bottom = 0.085 if colorkey else 0.0
    available = 1.0 - band_top - band_bottom

    structure_w = (xmax - xmin) * (1.0 + 2.0 * margin)
    structure_h = (ymax - ymin) * (1.0 + 2.0 * margin)
    view_h = max(structure_w / frame_aspect, structure_h / available)
    view_w = view_h * frame_aspect

    cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    ylo = cy - (band_bottom + 0.5 * available) * view_h
    ax.set_xlim(cx - 0.5 * view_w, cx + 0.5 * view_w)
    ax.set_ylim(ylo, ylo + view_h)
    ax.set_aspect("equal", adjustable="box")

    # --- text -------------------------------------------------------------
    # Measured rather than guessed: a hook is only as large as fits, so passing
    # a longer line shrinks the type instead of running it off the edge.
    scale = width_px / 1600.0
    bottom = 0.955
    if hook:
        artist = fig.text(0.045, bottom, textwrap.fill(hook, 30), ha="left",
                          va="top", fontsize=34 * scale, fontweight="bold",
                          color=COVER_INK, linespacing=1.15)
        bottom -= _fit_text(fig, artist, 0.90) + 0.025
    if subhook:
        artist = fig.text(0.045, bottom, textwrap.fill(subhook, 52), ha="left",
                          va="top", fontsize=15 * scale, color=COVER_INK_MUTED,
                          linespacing=1.3)
        _fit_text(fig, artist, 0.80)

    # --- minimal colour key ----------------------------------------------
    if colorkey:
        bar = fig.add_axes([0.045, 0.040, 0.20, 0.014])
        gradient = np.linspace(0.0, 1.0, 256).reshape(1, -1)
        bar.imshow(gradient, aspect="auto", cmap=st.FIELD_CMAP)
        bar.set_xticks([])
        bar.set_yticks([])
        for spine in bar.spines.values():
            spine.set_visible(False)
        bar.text(-0.04, 0.5, f"{vmin:.0f}", transform=bar.transAxes,
                 ha="right", va="center", fontsize=10 * scale,
                 color=COVER_INK_MUTED)
        bar.text(1.04, 0.5, f"{vmax:.0f} MPa", transform=bar.transAxes,
                 ha="left", va="center", fontsize=10 * scale,
                 color=COVER_INK_MUTED)
        bar.text(0.0, 2.3, "von Mises stress in the material",
                 transform=bar.transAxes, ha="left", va="bottom",
                 fontsize=10 * scale, color=COVER_INK_MUTED)

    fig.savefig(path, dpi=dpi, facecolor=COVER_SURFACE)
    plt.close(fig)
    print(f"  colour window {vmin:.1f} to {vmax:.1f} MPa "
          f"(p{low_percentile:g}-p{high_percentile:g} of the retained material; "
          f"true range {shown.min():.1f} to {shown.max():.1f} MPa)")
    return path


def auto_hook(case_dir: str) -> tuple:
    """A hook and sub-line the run's own summary supports without caveats."""
    from sparlab_viz.loaders import load_case as _load

    summary = _load(case_dir).summary
    fraction = summary.get("optimization_result", {}).get("volume_fraction")
    removed = f"{1.0 - fraction:.0%}" if fraction else "most"
    return (f"I deleted {removed} of this bracket.",
            "What is left is where the load actually goes.")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--case", default="results/aerospace_bracket")
    parser.add_argument("--load-case", default="down_limit")
    parser.add_argument("--output", default="docs/figures/cover_aerospace_bracket.png")
    parser.add_argument("--hook", default="What material does this structure "
                                          "actually need?")
    parser.add_argument("--subhook", default="")
    parser.add_argument("--hook-auto", action="store_true",
                        help="use a hook derived from the run's volume fraction")
    parser.add_argument("--no-hook", action="store_true")
    parser.add_argument("--no-colorkey", action="store_true")
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--pixels", default="1600x1200",
                        help="frame size, e.g. 1600x1200, 1200x1200, 1600x900")
    parser.add_argument("--low-percentile", type=float, default=2.0)
    parser.add_argument("--high-percentile", type=float, default=99.0)
    args = parser.parse_args(argv)

    hook, subhook = args.hook, args.subhook
    if args.hook_auto:
        hook, subhook = auto_hook(args.case)
    if args.no_hook:
        hook, subhook = "", ""

    try:
        w, h = (int(v) for v in args.pixels.lower().split("x"))
    except ValueError:
        print(f"--pixels must look like 1600x1200, got {args.pixels!r}",
              file=sys.stderr)
        return 1

    try:
        print(f"  wrote {build(args.case, args.load_case, args.output, hook, subhook, args.threshold, (w, h), args.low_percentile, args.high_percentile, not args.no_colorkey)}")
    except (ResultError, FileNotFoundError, KeyError) as error:
        print(f"could not build the cover: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
