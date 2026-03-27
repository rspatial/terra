"""
Pretty-print representations for SpatRaster, SpatVector, SpatExtent.

Mirrors R terra's ``show()`` / ``print()`` S4 methods (show.R).

Each ``repr_*`` function returns a multi-line string; the companion
``show_*`` functions print it.  ``register_reprs()`` monkey-patches
``__repr__`` onto the pybind11 types so that the REPL displays them
nicely automatically.
"""

from __future__ import annotations

import math
import os
from typing import Any, List, Optional

__all__ = [
    "repr_extent",
    "repr_raster",
    "repr_vector",
    "show",
    "register_reprs",
]

# ── Utilities ────────────────────────────────────────────────────────────────

_MAX_LYRS = 6          # maximum layers shown in name/min/max lines
_MAX_NAMES_W = 60      # total character budget for name columns
_MAX_SRC = 3           # maximum sources shown
_MAX_COLS_DF = 10      # max attribute columns shown in vector preview
_MAX_ROWS_DF = 3       # max attribute rows shown in vector preview


def _basename(path: str, n: int = 150) -> str:
    b = os.path.basename(path)
    if len(b) > n:
        b = b[:n] + "~"
    return b


def _crs_label(obj: Any) -> str:
    """Extract a short CRS label from WKT, falling back to PROJ.4."""
    try:
        wkt = obj.get_crs("wkt")
        if wkt:
            # First quoted string in WKT is usually the CRS name
            parts = wkt.split('"')
            if len(parts) >= 2:
                name = parts[1].strip()
                if name and name.lower() not in ("unknown", "unnamed", ""):
                    proj4 = obj.get_crs("proj4")
                    if proj4.startswith("+proj=longlat"):
                        return f"lon/lat {name}"
                    return name
    except Exception:
        pass
    try:
        p = obj.get_crs("proj4")
        return p if p else ""
    except Exception:
        return ""


def _fmt_num(v: float) -> str:
    if math.isnan(v) or math.isinf(v):
        return " ? "
    # Use up to 7 significant digits, strip trailing zeros
    s = f"{v:.7g}"
    return s


def _truncname(name: str, maxw: int) -> str:
    if len(name) <= maxw:
        return name
    mid = maxw // 2
    return name[:mid] + "~" + name[len(name) - mid + 1:]


# ── SpatExtent ──────────────────────────────────────────────────────────────

def repr_extent(e: Any) -> str:
    v = e.vector
    return (
        f"SpatExtent : {v[0]}, {v[1]}, {v[2]}, {v[3]}"
        "  (xmin, xmax, ymin, ymax)"
    )


# ── SpatRaster ───────────────────────────────────────────────────────────────

def repr_raster(r: Any) -> str:  # noqa: C901  (complex but mirrors R closely)
    lines: List[str] = []

    lines.append(f"class       : SpatRaster")

    nr, nc, nl = r.nrow(), r.ncol(), r.nlyr()
    lines.append(f"size        : {nr}, {nc}, {nl}  (nrow, ncol, nlyr)")

    rx, ry = r.res()
    lines.append(f"resolution  : {_fmt_num(rx)}, {_fmt_num(ry)}  (x, y)")

    e = r.extent.vector
    lines.append(
        f"extent      : {_fmt_num(e[0])}, {_fmt_num(e[1])}, "
        f"{_fmt_num(e[2])}, {_fmt_num(e[3])}  (xmin, xmax, ymin, ymax)"
    )

    lines.append(f"coord. ref. : {_crs_label(r)}")

    # ── Sources ──
    if r.hasValues:
        mem = r.inMemory
        fnames = r.filenames()
        nsrc = r.nsrc()

        if isinstance(mem, bool):
            mem = [mem] * nsrc
        if isinstance(mem, (list, tuple)) and len(mem) < nsrc:
            mem = list(mem) + [False] * (nsrc - len(mem))

        srcs: List[str] = []
        for i, f in enumerate(fnames):
            if isinstance(mem, (list, tuple)) and i < len(mem) and mem[i]:
                srcs.append("memory")
            else:
                srcs.append(_basename(f) if f else "memory")

        if all(s == "memory" for s in srcs):
            lines.append("source(s)   : memory")
        elif nsrc == 1:
            lines.append(f"source      : {srcs[0]}")
        else:
            for i, s in enumerate(srcs[:_MAX_SRC]):
                label = "sources     :" if i == 0 else "             "
                lines.append(f"{label} {s}")
            if nsrc > _MAX_SRC:
                rem = nsrc - _MAX_SRC
                if rem == 1:
                    lines.append(f"              {srcs[_MAX_SRC]}")
                else:
                    lines.append(f"              ... and {rem} more sources")

        # ── Names + min/max ──
        raw_names = list(r.names)
        disp_names = raw_names[:_MAX_LYRS]
        if nl > _MAX_LYRS:
            disp_names.append("...")

        maxw = max(1, _MAX_NAMES_W // max(1, min(_MAX_LYRS, len(disp_names))))
        disp_names = [_truncname(n, maxw) for n in disp_names]

        has_range = r.hasRange
        if isinstance(has_range, bool):
            has_range = [has_range] * nl
        has_any = any(has_range[:_MAX_LYRS])

        if has_any:
            rmins = list(r.range_min)
            rmaxs = list(r.range_max)

            def _fmt_val(v: float, has: bool, suffix: str = "") -> str:
                if not has:
                    return " ? "
                return suffix + _fmt_num(v)

            minv = [_fmt_val(rmins[i], has_range[i]) for i in range(min(nl, _MAX_LYRS))]
            maxv = [_fmt_val(rmaxs[i], has_range[i]) for i in range(min(nl, _MAX_LYRS))]
            if nl > _MAX_LYRS:
                minv.append("...")
                maxv.append("...")

            # column widths
            w = [
                max(len(n), len(mn), len(mx))
                for n, mn, mx in zip(disp_names, minv, maxv)
            ]
            row_n = "  ".join(n.rjust(w[i]) for i, n in enumerate(disp_names))
            row_min = "  ".join(v.rjust(w[i]) for i, v in enumerate(minv))
            row_max = "  ".join(v.rjust(w[i]) for i, v in enumerate(maxv))

            lbl = "name        :" if nl == 1 else "names       :"
            lines.append(f"{lbl} {row_n}")
            lbl = "min value   :" if nl == 1 else "min values  :"
            lines.append(f"{lbl} {row_min}")
            lbl = "max value   :" if nl == 1 else "max values  :"
            lines.append(f"{lbl} {row_max}")
        else:
            lbl = "name        :" if nl == 1 else "names       :"
            lines.append(f"{lbl} {',  '.join(disp_names)}")

        # ── Units ──
        try:
            units = r.getUnit()
            unique_units = list(dict.fromkeys(u for u in units if u))
            if unique_units:
                if len(unique_units) == 1:
                    lines.append(f"unit        : {unique_units[0]}")
                else:
                    uts = units[:_MAX_LYRS]
                    if nl > _MAX_LYRS:
                        uts.append("...")
                    lines.append(f"unit        : {',  '.join(uts)}")
        except Exception:
            pass

    # ── Time ──
    if r.hasTime:
        try:
            tms = r.time
            step = r.timestep
            tz = r.timezone

            step_labels = {
                "yearmonths": "time (ymnts)",
                "months":     "time (mnts) ",
                "years":      "time (years)",
                "days":       "time (days) ",
                "raw":        "time (raw)  ",
            }
            label = step_labels.get(step, "time        ")

            if tms:
                tmin, tmax = min(tms), max(tms)
                n_steps = len(set(tms))
                if tmin == tmax:
                    t_str = str(tmin)
                else:
                    t_str = f"{tmin} to {tmax} ({n_steps} steps)"
                if step == "seconds" and tz:
                    t_str += f" {tz}"
                lines.append(f"{label}: {t_str}")
        except Exception:
            pass

    # ── Depth ──
    if r.hasDepth:
        try:
            depths = r.getDepth()
            dname = r.getDepthName()
            dunit = r.getDepthUnit()

            if dname == "depth":
                prefix = f"[{dunit}]: " if dunit else ""
            else:
                if dunit and dunit != "unknown":
                    prefix = f"{dname} [{dunit}]: "
                else:
                    prefix = f"{dname}: "

            if depths:
                udepths = sorted(set(depths))
                if len(udepths) > 1:
                    depth_str = (
                        f"{udepths[0]} to {udepths[-1]} "
                        f"({prefix}{len(udepths)} steps)"
                    )
                else:
                    depth_str = f"{prefix}{udepths[0]}"
                lines.append(f"depth       : {depth_str}")
        except Exception:
            pass

    return "\n".join(lines)


# ── SpatVector ───────────────────────────────────────────────────────────────

def _fmt_df_preview(df_dict: dict, nrow: int) -> str:
    """
    Render a simple attribute-table preview (first _MAX_ROWS_DF rows,
    first _MAX_COLS_DF columns) as a plain-text block, like R's printDF.
    """
    if not df_dict:
        return ""

    cols = list(df_dict.keys())[:_MAX_COLS_DF]
    extra_cols = len(df_dict) - len(cols)
    rows_shown = min(nrow, _MAX_ROWS_DF)

    if rows_shown == 0:
        return ""

    # Build type row
    def _pytype(col: list) -> str:
        if not col:
            return "<?>"
        v = next((x for x in col if x is not None), None)
        if v is None:
            return "<NA>"
        if isinstance(v, bool):
            return "<bool>"
        if isinstance(v, int):
            return "<int>"
        if isinstance(v, float):
            return "<num>"
        return "<chr>"

    type_row = [_pytype(df_dict[c]) for c in cols]

    # Build value rows
    value_rows = []
    for i in range(rows_shown):
        row = []
        for c in cols:
            vals = df_dict[c]
            v = vals[i] if i < len(vals) else None
            if v is None:
                row.append("NA")
            elif isinstance(v, float):
                row.append(f"{v:.6g}")
            else:
                row.append(str(v))
        value_rows.append(row)

    all_rows = [type_row] + value_rows
    if nrow > rows_shown:
        all_rows.append(["..."] * len(cols))

    # Column widths
    widths = [
        max(len(c), max(len(r[i]) for r in all_rows))
        for i, c in enumerate(cols)
    ]

    def _fmt_row(tag: str, row: list) -> str:
        cells = "  ".join(v.rjust(widths[i]) for i, v in enumerate(row))
        return f"{tag} {cells}"

    out_lines: List[str] = []
    out_lines.append(_fmt_row("type        :", type_row))
    if value_rows:
        out_lines.append(_fmt_row("values      :", value_rows[0]))
        for vrow in value_rows[1:]:
            out_lines.append(_fmt_row("             ", vrow))
    if nrow > rows_shown:
        out_lines.append("             ...")

    if extra_cols:
        out_lines.append(f"              (and {extra_cols} more attribute{'s' if extra_cols>1 else ''})")

    # Column name header above type row
    header_tag = "names       :"
    hcells = "  ".join(c.rjust(widths[i]) for i, c in enumerate(cols))
    out_lines.insert(0, f"{header_tag} {hcells}")

    return "\n".join(out_lines)


def repr_vector(v: Any) -> str:
    lines: List[str] = []

    lines.append(" class       : SpatVector")

    geom = v.type()
    lines.append(f" geometry    : {geom}")

    nr, nc = v.nrow(), v.ncol()
    lines.append(f" dimensions  : {nr}, {nc}  (geometries, attributes)")

    try:
        e = v.extent().vector
        lines.append(
            f" extent      : {_fmt_num(e[0])}, {_fmt_num(e[1])}, "
            f"{_fmt_num(e[2])}, {_fmt_num(e[3])}  (xmin, xmax, ymin, ymax)"
        )
    except Exception:
        pass

    src = getattr(v, "source", "")
    lyr = getattr(v, "layer", "")
    if src:
        b = _basename(src)
        sans_ext = os.path.splitext(b)[0]
        if lyr and lyr != sans_ext:
            lines.append(f" source      : {b} ({lyr})")
        else:
            lines.append(f" source      : {b}")

    lines.append(f" coord. ref. : {_crs_label(v)}")

    if nc > 0 and nr > 0:
        try:
            df = v.getDF()
            preview = _fmt_df_preview(df, nr)
            if preview:
                lines.append(preview)
        except Exception:
            pass

    return "\n".join(lines)


# ── Generic dispatcher ────────────────────────────────────────────────────────

def show(obj: Any) -> None:
    """Print a repr for any supported terra object."""
    from ._terra import SpatExtent, SpatRaster, SpatVector

    if isinstance(obj, SpatRaster):
        print(repr_raster(obj))
    elif isinstance(obj, SpatVector):
        print(repr_vector(obj))
    elif isinstance(obj, SpatExtent):
        print(repr_extent(obj))
    else:
        print(repr(obj))


# ── Monkey-patch __repr__ onto pybind11 types ────────────────────────────────

def register_reprs() -> None:
    """
    Attach ``__repr__`` and ``__str__`` to the C++ pybind11 classes so that
    the REPL shows nice output automatically.

    Called once at import time from ``__init__.py``.
    """
    from ._terra import SpatExtent, SpatRaster, SpatVector

    SpatRaster.__repr__ = lambda self: repr_raster(self)   # type: ignore[attr-defined]
    SpatRaster.__str__  = lambda self: repr_raster(self)   # type: ignore[attr-defined]
    SpatVector.__repr__ = lambda self: repr_vector(self)   # type: ignore[attr-defined]
    SpatVector.__str__  = lambda self: repr_vector(self)   # type: ignore[attr-defined]
    SpatExtent.__repr__ = lambda self: repr_extent(self)   # type: ignore[attr-defined]
    SpatExtent.__str__  = lambda self: repr_extent(self)   # type: ignore[attr-defined]
