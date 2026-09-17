"""SparLab result loading and plotting.

This package reads only what the C++ apps wrote: `mesh.json`, `summary.json`
and the CSV tables in a result directory. It performs no analysis of its own, so
every number in a figure traces back to a solver run. The one exception is
clearly-labelled derived quantities (mass from volume and density, error ratios),
and those are computed from values in the same summary.

Layout
------
loaders   reading mesh.json / summary.json / the CSV tables
style     one place for figure size, fonts, colour maps and the save helper
fields    turning element and nodal fields into matplotlib collections
plots     the individual figures the documentation uses
"""

from .loaders import (
    CaseResults,
    load_case,
    load_csv,
    load_json,
    load_mesh,
    load_summary,
)
from .style import (
    DENSITY_CMAP,
    FIELD_CMAP,
    SIGNED_CMAP,
    apply_style,
    figure,
    save_figure,
)
from .fields import (
    element_collection,
    nodal_tripcolor,
    deformed_collection,
    add_colorbar,
    set_domain_limits,
)

__all__ = [
    "CaseResults",
    "load_case",
    "load_csv",
    "load_json",
    "load_mesh",
    "load_summary",
    "DENSITY_CMAP",
    "FIELD_CMAP",
    "SIGNED_CMAP",
    "apply_style",
    "figure",
    "save_figure",
    "element_collection",
    "nodal_tripcolor",
    "deformed_collection",
    "add_colorbar",
    "set_domain_limits",
]

__version__ = "1.0.0"
