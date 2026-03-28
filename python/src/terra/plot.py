"""
Plotting for SpatRaster objects.

Provides :func:`plot` (single or multi-layer raster) and :func:`plot_rgb`
(composite colour image) using **matplotlib** as the rendering backend.

Quick start::

    import terra as pt
    import matplotlib.pyplot as plt

    r = pt.rast("elevation.tif")
    pt.plot(r)
    plt.show()
"""

from __future__ import annotations

import math
import warnings
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np

from ._terra import SpatOptions, SpatRaster

__all__ = ["plot", "plot_rgb"]

# ── Default palette ──────────────────────────────────────────────────────────

def _default_palette(n: int = 255) -> List[str]:
    """
    Return the default terra colour ramp (a reverse of matplotlib's terrain_r).

    Args:
        n: Number of colours.

    Returns:
        List of hex colour strings.
    """
    import matplotlib.cm as cm
    import matplotlib.colors as mc
    cmap = cm.get_cmap("terrain_r", n)
    return [mc.to_hex(cmap(i)) for i in range(n)]


# ── Value extraction ──────────────────────────────────────────────────────────

def _get_layer_array(r: SpatRaster, lyr: int) -> np.ndarray:
    """
    Read one layer (0-based index) as a 2-D float64 array (rows × cols).

    NA / Inf values become ``np.nan``.
    """
    opt = SpatOptions()
    flat = r.getValues(lyr, opt)
    arr = np.array(flat, dtype=np.float64).reshape(r.nrow(), r.ncol())
    arr[~np.isfinite(arr)] = np.nan
    return arr


def _downsample_array(arr: np.ndarray, max_cells: int) -> np.ndarray:
    """
    Thin a 2-D array by a uniform stride so that ncell ≤ max_cells.

    Args:
        arr: 2-D numpy array (nrow × ncol).
        max_cells: Target maximum number of cells.

    Returns:
        Thinned array.
    """
    nr, nc = arr.shape
    total = nr * nc
    if total <= max_cells:
        return arr
    factor = math.ceil(math.sqrt(total / max_cells))
    return arr[::factor, ::factor]


# ── Colour helpers ────────────────────────────────────────────────────────────

def _hex_to_rgba(hex_colors: Sequence[str]) -> np.ndarray:
    """Convert a sequence of hex colour strings to an (N, 4) float array."""
    import matplotlib.colors as mc
    return np.array([mc.to_rgba(c) for c in hex_colors], dtype=np.float64)


def _breaks_equal_interval(
    values: np.ndarray, n_bins: int, range_vals: Tuple[float, float]
) -> np.ndarray:
    """Compute ``n_bins + 1`` equally-spaced break points over *range_vals*."""
    lo, hi = range_vals
    if lo == hi:
        return np.array([lo - 0.5, lo + 0.5])
    return np.linspace(lo, hi, n_bins + 1)


def _auto_digits(value_range: float) -> int:
    """Return appropriate decimal digit count for a legend given *value_range*."""
    if value_range == 0 or not math.isfinite(value_range):
        return 0
    return max(0, -math.floor(math.log10(value_range / 10)))


# ── Continuous pipeline ───────────────────────────────────────────────────────

def _continuous_image(
    arr: np.ndarray,
    palette: List[str],
    range_vals: Optional[Tuple[Optional[float], Optional[float]]] = None,
    fill_range: bool = False,
) -> Tuple[np.ndarray, Tuple[float, float], int]:
    """
    Map a continuous 2-D array to an RGBA image using *palette*.

    Args:
        arr: 2-D float array (rows × cols), NaN for no-data.
        palette: Ordered list of hex colour strings.
        range_vals: ``(vmin, vmax)`` clipping range.  ``None`` elements are
            filled from the data.
        fill_range: If True, clamp out-of-range values to the range endpoints
            rather than masking them.

    Returns:
        ``(rgba_image, (vmin, vmax), n_digits)`` where *rgba_image* has shape
        ``(rows, cols, 4)``.
    """
    import matplotlib.colors as mc

    valid = arr[np.isfinite(arr)]
    if valid.size == 0:
        rgba = np.zeros((*arr.shape, 4), dtype=np.float64)
        return rgba, (np.nan, np.nan), 0

    # Resolve data range
    data_min, data_max = float(np.nanmin(valid)), float(np.nanmax(valid))
    vmin: float = data_min
    vmax: float = data_max

    if range_vals is not None:
        lo, hi = range_vals
        if fill_range:
            if lo is not None and np.isfinite(lo):
                arr = np.where(arr < lo, lo, arr)
                vmin = lo
            else:
                vmin = data_min
            if hi is not None and np.isfinite(hi):
                arr = np.where(arr > hi, hi, arr)
                vmax = hi
            else:
                vmax = data_max
        else:
            vmin = lo if (lo is not None and np.isfinite(lo)) else data_min
            vmax = hi if (hi is not None and np.isfinite(hi)) else data_max

    n = len(palette)
    cmap = mc.ListedColormap(palette)

    if vmin == vmax:
        # Only one unique value — use the middle colour
        norm = mc.Normalize(vmin=vmin - 0.5, vmax=vmax + 0.5)
    else:
        norm = mc.Normalize(vmin=vmin, vmax=vmax)

    rgba = cmap(norm(np.ma.masked_invalid(arr)))
    # Restore transparency where data is NaN
    nan_mask = ~np.isfinite(arr)
    rgba[nan_mask, 3] = 0.0

    digits = _auto_digits(vmax - vmin)
    return rgba, (vmin, vmax), digits


# ── Classes pipeline ──────────────────────────────────────────────────────────

def _classes_image(
    arr: np.ndarray,
    palette: List[str],
    levels: Optional[List[float]] = None,
) -> Tuple[np.ndarray, List[float], List[str]]:
    """
    Map discrete class values in *arr* to colours from *palette*.

    Args:
        arr: 2-D float array.  Values not listed in *levels* become
            transparent.
        palette: List of hex colour strings.
        levels: Ordered list of numeric values to map.  If None, the unique
            non-NaN values in *arr* are used.

    Returns:
        ``(rgba_image, levels, hex_fill_colours)``
    """
    import matplotlib.colors as mc

    valid = arr[np.isfinite(arr)]
    if valid.size == 0:
        rgba = np.zeros((*arr.shape, 4), dtype=np.float64)
        return rgba, [], []

    if levels is None:
        levels = sorted(set(float(v) for v in valid))

    n_lvl = len(levels)
    n_col = len(palette)
    if n_lvl == 1:
        cols = [palette[-1]]
    elif n_lvl <= n_col:
        idx = np.round(np.linspace(0, n_col - 1, n_lvl)).astype(int)
        cols = [palette[i] for i in idx]
    else:
        cols = [palette[i % n_col] for i in range(n_lvl)]

    level_to_col = {lv: mc.to_rgba(c) for lv, c in zip(levels, cols)}

    rgba = np.zeros((*arr.shape, 4), dtype=np.float64)
    for lv, color in level_to_col.items():
        mask = arr == lv
        rgba[mask] = color

    return rgba, levels, cols


# ── Factor pipeline ───────────────────────────────────────────────────────────

def _cats_to_dict(cats_obj: Any) -> Tuple[List[float], List[str]]:
    """
    Extract ``(ids, labels)`` from a SpatCategories object.

    Returns a list of numeric IDs and a parallel list of label strings.
    """
    df = cats_obj.d
    # SpatDataFrame columns come back as a dict-like; try common key patterns
    try:
        ids = list(df.get_column(0))
        label_col = cats_obj.index if hasattr(cats_obj, "index") else 1
        labels = [str(v) for v in df.get_column(label_col)]
    except Exception:
        ids, labels = [], []
    return ids, labels


def _coltab_to_hex(coltab_obj: Any) -> Tuple[List[float], List[str]]:
    """
    Convert a SpatDataFrame colour table to ``(ids, hex_colours)``.

    The SpatDataFrame has columns: id, R, G, B, A (values 0–255).
    """
    df = coltab_obj
    try:
        ids = list(df.get_column(0))
        r = np.array(df.get_column(1), dtype=np.uint8)
        g = np.array(df.get_column(2), dtype=np.uint8)
        b = np.array(df.get_column(3), dtype=np.uint8)
        a = np.array(df.get_column(4), dtype=np.uint8)
        hexcols = [
            "#{:02x}{:02x}{:02x}{:02x}".format(ri, gi, bi, ai)
            for ri, gi, bi, ai in zip(r, g, b, a)
        ]
    except Exception:
        ids, hexcols = [], []
    return ids, hexcols


def _factor_image(
    arr: np.ndarray,
    cats_obj: Any,
    coltab_obj: Optional[Any],
    palette: List[str],
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Map a factor-coded raster to RGBA colours, using the category table.

    Args:
        arr: 2-D integer-valued float array.
        cats_obj: SpatCategories for this layer.
        coltab_obj: Optional colour-table SpatDataFrame for this layer.
        palette: Fallback colour palette.

    Returns:
        ``(rgba_image, legend_labels, hex_fill_colours)``
    """
    import matplotlib.colors as mc

    ids, labels = _cats_to_dict(cats_obj)
    if not ids:
        return np.zeros((*arr.shape, 4)), [], []

    if coltab_obj is not None:
        ct_ids, ct_hex = _coltab_to_hex(coltab_obj)
        id_to_color = {ct_ids[i]: ct_hex[i] for i in range(len(ct_ids))}
    else:
        n_lvl = len(ids)
        n_col = len(palette)
        if n_lvl <= n_col:
            idxs = np.round(np.linspace(0, n_col - 1, n_lvl)).astype(int)
            selected = [palette[i] for i in idxs]
        else:
            selected = [palette[i % n_col] for i in range(n_lvl)]
        id_to_color = {ids[i]: selected[i] for i in range(len(ids))}

    fill_cols = [id_to_color.get(iid, "#ffffff00") for iid in ids]

    rgba = np.zeros((*arr.shape, 4), dtype=np.float64)
    for iid, color in id_to_color.items():
        mask = arr == float(iid)
        rgba[mask] = mc.to_rgba(color)

    return rgba, labels, fill_cols


# ── Interval pipeline ─────────────────────────────────────────────────────────

def _interval_image(
    arr: np.ndarray,
    breaks: Union[np.ndarray, List[float]],
    palette: List[str],
) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Map values to colours according to user-defined break intervals.

    Each colour in *palette* is assigned to the half-open interval
    ``[breaks[i], breaks[i+1])``, matching R's ``cut(..., right=FALSE)``.

    Args:
        arr: 2-D float array.
        breaks: Monotonically increasing break points.  ``len(breaks)`` must
            equal ``len(palette) + 1``.
        palette: One colour per interval.

    Returns:
        ``(rgba_image, interval_labels, hex_fill_colours)``
    """
    import matplotlib.colors as mc

    breaks = np.asarray(breaks, dtype=np.float64)
    n_bins = len(breaks) - 1
    cols = palette[:n_bins] if len(palette) >= n_bins else palette + [palette[-1]] * (n_bins - len(palette))

    bin_idx = np.digitize(arr, breaks[:-1]) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    rgba = np.zeros((*arr.shape, 4), dtype=np.float64)
    for i, c in enumerate(cols):
        mask = (bin_idx == i) & np.isfinite(arr)
        rgba[mask] = mc.to_rgba(c)

    labels = [f"{breaks[i]:.4g} – {breaks[i+1]:.4g}" for i in range(n_bins)]
    return rgba, labels, cols


# ── Colortable pipeline ───────────────────────────────────────────────────────

def _colortable_image(arr: np.ndarray, coltab_obj: Any) -> np.ndarray:
    """
    Map integer cell values to RGBA colours using an embedded colour table.

    Args:
        arr: 2-D integer-valued float array.
        coltab_obj: SpatDataFrame colour table (id, R, G, B, A columns,
            values 0–255).

    Returns:
        RGBA image array (rows × cols × 4).
    """
    import matplotlib.colors as mc

    ct_ids, ct_hex = _coltab_to_hex(coltab_obj)
    if not ct_ids:
        return np.zeros((*arr.shape, 4))

    id_to_color = {float(iid): mc.to_rgba(c) for iid, c in zip(ct_ids, ct_hex)}
    rgba = np.zeros((*arr.shape, 4), dtype=np.float64)
    for iid, color in id_to_color.items():
        mask = arr == iid
        rgba[mask] = color
    return rgba


# ── RGB pipeline ──────────────────────────────────────────────────────────────

def _stretch_band(band: np.ndarray, method: Optional[str]) -> np.ndarray:
    """
    Stretch band values to the [0, 1] range.

    Args:
        band: 2-D float array.
        method: ``"lin"`` for linear 2%–98% percentile stretch, ``"hist"``
            for histogram equalisation, or ``None`` for simple min-max.

    Returns:
        Array with values in [0, 1].
    """
    valid = band[np.isfinite(band)]
    if valid.size == 0:
        return np.zeros_like(band)
    if method in ("lin", "linear"):
        lo, hi = np.percentile(valid, [2, 98])
    else:
        lo, hi = float(np.nanmin(valid)), float(np.nanmax(valid))
    if lo == hi:
        return np.zeros_like(band)
    out = np.clip((band - lo) / (hi - lo), 0.0, 1.0)
    out[~np.isfinite(band)] = 0.0
    return out


def _rgb_image(
    r: SpatRaster,
    rgb_bands: Sequence[int],
    scale: float = 255.0,
    stretch: Optional[str] = None,
    na_color: str = "white",
) -> np.ndarray:
    """
    Compose an RGBA image from three (or four) raster bands.

    Args:
        r: SpatRaster with at least three layers.
        rgb_bands: 0-based layer indices for Red, Green, Blue (and optionally
            Alpha).
        scale: Maximum value used for normalisation (typically 255 for 8-bit
            data, or 65535 for 16-bit data).
        stretch: Optional contrast stretch.  ``"lin"`` applies a linear
            2%–98% percentile stretch; ``"hist"`` applies histogram
            equalisation.
        na_color: Hex or named colour used for cells where any band is NA.

    Returns:
        RGBA array (rows × cols × 4) with values in [0, 1].
    """
    import matplotlib.colors as mc

    bands = []
    for b in rgb_bands[:3]:
        arr = _get_layer_array(r, b)
        if stretch:
            arr = _stretch_band(arr, stretch)
        else:
            arr = np.clip(arr / scale, 0.0, 1.0)
            arr[~np.isfinite(arr)] = 0.0
        bands.append(arr)

    nr, nc = bands[0].shape
    rgba = np.ones((nr, nc, 4), dtype=np.float64)
    rgba[:, :, 0] = bands[0]
    rgba[:, :, 1] = bands[1]
    rgba[:, :, 2] = bands[2]

    if len(rgb_bands) >= 4:
        alpha_arr = _get_layer_array(r, rgb_bands[3])
        rgba[:, :, 3] = np.clip(alpha_arr / scale, 0.0, 1.0)
        rgba[:, :, 3][~np.isfinite(alpha_arr)] = 0.0

    # Mask cells where any band was NA
    na_mask = np.zeros((nr, nc), dtype=bool)
    for b in rgb_bands[:3]:
        arr = _get_layer_array(r, b)
        na_mask |= ~np.isfinite(arr)
    na_rgba = mc.to_rgba(na_color)
    rgba[na_mask] = na_rgba

    return rgba


# ── Legend helpers ────────────────────────────────────────────────────────────

def _add_continuous_legend(
    ax: Any,
    cmap: Any,
    norm: Any,
    n_ticks: int = 5,
    digits: int = 2,
    title: str = "",
) -> None:
    """Add a colourbar (continuous legend) to *ax*."""
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm

    fig = ax.get_figure()
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(title)
    ticks = np.linspace(norm.vmin, norm.vmax, n_ticks)
    cb.set_ticks(ticks)
    fmt = f"{{:.{digits}f}}"
    cb.set_ticklabels([fmt.format(t) for t in ticks])


def _add_class_legend(
    ax: Any,
    labels: List[str],
    colors: List[str],
    title: str = "",
    reverse: bool = False,
) -> None:
    """Add a categorical legend (patches) to *ax*."""
    import matplotlib.patches as mpatches

    if reverse:
        labels = labels[::-1]
        colors = colors[::-1]
    patches = [
        mpatches.Patch(facecolor=c, label=l)
        for c, l in zip(colors, labels)
    ]
    ax.legend(
        handles=patches,
        title=title,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        borderaxespad=0,
        frameon=True,
    )


# ── Axes / title helpers ──────────────────────────────────────────────────────

def _setup_axes(
    ax: Any,
    ext: Sequence[float],
    axes: bool,
    lonlat: bool,
) -> None:
    """Configure axis ticks and labels for a spatial plot."""
    xmin, xmax, ymin, ymax = ext
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    if not axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis("off")
    else:
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.tick_params(labelsize=8)


def _get_nrnc(nr: Optional[int], nc: Optional[int], nl: int) -> Tuple[int, int]:
    """Determine the subplot grid dimensions for *nl* layers."""
    if nc is not None and nr is not None:
        return int(nr), int(nc)
    if nc is not None:
        return math.ceil(nl / int(nc)), int(nc)
    if nr is not None:
        return int(nr), math.ceil(nl / int(nr))
    nc_auto = math.ceil(math.sqrt(nl))
    nr_auto = math.ceil(nl / nc_auto)
    return int(nr_auto), int(nc_auto)


# ── Infer plot type ───────────────────────────────────────────────────────────

def _infer_type(r: SpatRaster, lyr: int) -> str:
    """
    Decide the plot type for layer *lyr* (0-based), mirroring R's auto-detection.

    Returns one of ``"continuous"``, ``"factor"``, ``"colortable"``,
    ``"depends"``, or ``"rgb"``.
    """
    rgb = r.getRGB()
    if rgb and lyr in rgb:
        return "rgb"
    vt = r.valueType(False)
    has_cats = r.hasCategories()
    has_col = r.hasColors()
    if has_cats[lyr] and has_col[lyr]:
        return "factor"
    if has_col[lyr]:
        return "colortable"
    if has_cats[lyr]:
        return "factor"
    if vt[lyr] == 3:   # boolean
        return "factor"
    return "depends"


# ── Single-layer plot ─────────────────────────────────────────────────────────

def _plot_one_layer(
    r: SpatRaster,
    lyr: int,
    ax: Any,
    palette: List[str],
    type: str,
    zlim: Optional[Tuple[Optional[float], Optional[float]]],
    clamp: bool,
    breaks: Optional[Union[Sequence[float], np.ndarray]],
    levels: Optional[List[Any]],
    legend: bool,
    na_color: Optional[str],
    axes: bool,
    title: str,
    max_cell: int,
    alpha: float,
    smooth: bool,
) -> None:
    """
    Render a single raster layer into *ax*.

    This is the workhorse called by :func:`plot` for each panel.
    """
    import matplotlib.colors as mc

    ext_v = r.extent.vector                   # (xmin, xmax, ymin, ymax)
    lonlat = r.isLonLat()

    arr = _get_layer_array(r, lyr)
    arr = _downsample_array(arr, max_cell)

    interp = "bilinear" if smooth else "nearest"

    # ── detect type ──────────────────────────────────────────────────────────
    if type == "depends":
        valid_u = np.unique(arr[np.isfinite(arr)])
        if len(valid_u) < 9:
            type = "classes"
        else:
            type = "continuous"

    # ── colour pipeline ──────────────────────────────────────────────────────
    if type in ("continuous", "interval") and breaks is not None:
        breaks_arr = np.asarray(breaks, dtype=np.float64)
        rgba, labels, fill_cols = _interval_image(arr, breaks_arr, palette)
        _setup_axes(ax, ext_v, axes, lonlat)
        ax.imshow(rgba, extent=ext_v, origin="upper",
                  interpolation=interp, aspect="auto", alpha=alpha)
        if legend:
            _add_class_legend(ax, labels, fill_cols, title=title)

    elif type == "continuous":
        rgba, (vmin, vmax), digits = _continuous_image(
            arr, palette, range_vals=zlim, fill_range=clamp
        )
        _setup_axes(ax, ext_v, axes, lonlat)
        n = len(palette)
        cmap = mc.ListedColormap(palette)
        norm = mc.Normalize(
            vmin=vmin if np.isfinite(vmin) else 0,
            vmax=vmax if np.isfinite(vmax) else 1,
        )
        if not np.isnan(vmin) and not np.isnan(vmax):
            ax.imshow(rgba, extent=ext_v, origin="upper",
                      interpolation=interp, aspect="auto", alpha=alpha)
            if legend:
                _add_continuous_legend(ax, cmap, norm, digits=digits, title=title)
        else:
            ax.imshow(rgba, extent=ext_v, origin="upper",
                      interpolation=interp, aspect="auto")

    elif type == "classes":
        lv = [float(v) for v in levels] if levels is not None else None
        rgba, lv_used, fill_cols = _classes_image(arr, palette, levels=lv)
        _setup_axes(ax, ext_v, axes, lonlat)
        ax.imshow(rgba, extent=ext_v, origin="upper",
                  interpolation=interp, aspect="auto", alpha=alpha)
        if legend:
            str_labels = [str(v) for v in lv_used]
            _add_class_legend(ax, str_labels, fill_cols, title=title)

    elif type == "factor":
        cats_list = r.getCategories()
        cats_obj = cats_list[lyr] if cats_list else None
        coltabs = r.getColors()
        coltab_obj = coltabs[lyr] if coltabs and r.hasColors()[lyr] else None
        if cats_obj is not None:
            rgba, leg_labels, fill_cols = _factor_image(
                arr, cats_obj, coltab_obj, palette
            )
        else:
            rgba, lv_used, fill_cols = _classes_image(arr, palette)
            leg_labels = [str(v) for v in lv_used]
        _setup_axes(ax, ext_v, axes, lonlat)
        ax.imshow(rgba, extent=ext_v, origin="upper",
                  interpolation=interp, aspect="auto", alpha=alpha)
        if legend:
            _add_class_legend(ax, leg_labels, fill_cols, title=title)

    elif type == "colortable":
        coltabs = r.getColors()
        coltab_obj = coltabs[lyr] if coltabs else None
        if coltab_obj is not None:
            rgba = _colortable_image(arr, coltab_obj)
        else:
            rgba, _, _ = _classes_image(arr, palette)
        _setup_axes(ax, ext_v, axes, lonlat)
        ax.imshow(rgba, extent=ext_v, origin="upper",
                  interpolation=interp, aspect="auto", alpha=alpha)

    elif type == "rgb":
        rgb_idxs = r.getRGB()
        if not rgb_idxs:
            rgb_idxs = [0, 1, 2]
        rgba = _rgb_image(r, rgb_idxs, stretch=None,
                          na_color=na_color or "white")
        _setup_axes(ax, ext_v, axes, lonlat)
        ax.imshow(rgba, extent=ext_v, origin="upper",
                  interpolation=interp, aspect="auto")

    # ── NA colour overlay ────────────────────────────────────────────────────
    if na_color is not None and type not in ("rgb", "colortable"):
        _add_na_overlay(ax, arr, na_color, ext_v, interp, alpha)

    # ── title ────────────────────────────────────────────────────────────────
    if title:
        ax.set_title(title, fontsize=9)


def _add_na_overlay(
    ax: Any,
    arr: np.ndarray,
    na_color: str,
    ext: Sequence[float],
    interp: str,
    alpha: float,
) -> None:
    """Draw a solid overlay colour for NA cells."""
    import matplotlib.colors as mc
    import numpy as np

    na_rgba = mc.to_rgba(na_color)
    overlay = np.zeros((*arr.shape, 4), dtype=np.float64)
    mask = ~np.isfinite(arr)
    overlay[mask] = na_rgba
    ax.imshow(overlay, extent=ext, origin="upper",
              interpolation=interp, aspect="auto")


# ── Public API ────────────────────────────────────────────────────────────────

def plot(
    r: SpatRaster,
    y: Union[int, List[int], str, List[str], None] = None,
    *,
    col: Optional[Union[List[str], str]] = None,
    type: Optional[str] = None,
    legend: bool = True,
    axes: bool = True,
    zlim: Optional[Tuple[Optional[float], Optional[float]]] = None,
    clamp: bool = False,
    levels: Optional[List[Any]] = None,
    breaks: Optional[Sequence[float]] = None,
    na_color: Optional[str] = "white",
    alpha: float = 1.0,
    smooth: bool = False,
    maxcell: int = 500_000,
    nc: Optional[int] = None,
    nr: Optional[int] = None,
    maxnl: int = 16,
    main: Optional[Union[str, List[str]]] = None,
    ax: Optional[Any] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs: Any,
) -> Any:
    """
    Plot one or more layers of a SpatRaster.

    When ``y`` selects a single layer the result is drawn in *ax* (or a newly
    created Axes).  When multiple layers are selected, a grid of subplots is
    created and all Axes are returned as a 2-D numpy array.

    Args:
        r: SpatRaster to plot.
        y: Layer selector.  Can be:

            * ``None`` — plot all layers (up to ``maxnl``).
            * ``int`` — 0-based layer index (Python convention; ``-1`` is last).
            * ``list[int]`` — multiple 0-based layer indices.
            * ``str`` or ``list[str]`` — layer name(s).

        col: Colour palette as a list of hex strings, or a matplotlib
            colormap name (e.g. ``"viridis"``).  Defaults to the terra
            terrain palette.
        type: Plot type.  One of ``"continuous"``, ``"classes"``,
            ``"factor"``, ``"interval"``, ``"colortable"``, ``"rgb"``.
            If None, the type is inferred from the data.
        legend: If True (default), draw a legend or colourbar.
        axes: If True (default), draw coordinate axis ticks and labels.
        zlim: ``(vmin, vmax)`` display range for continuous data.  ``None``
            elements are derived from the data.
        clamp: If True, values outside *zlim* are clamped to the range
            endpoints rather than shown as NA.
        levels: Explicit numeric values to use as class levels.
        breaks: Cut-point values for interval classification.  When supplied,
            *type* is automatically set to ``"interval"``.
        na_color: Colour for NA cells.  Use ``None`` to make them
            transparent.
        alpha: Overall opacity for the raster image (0–1).
        smooth: If True, apply bilinear interpolation when rendering
            (relevant only for display, not the data).
        maxcell: Maximum number of cells to render.  The raster is thinned
            by a uniform stride when ``ncell(r) > maxcell``.
        nc: Number of subplot columns (multi-layer only).
        nr: Number of subplot rows (multi-layer only).
        maxnl: Maximum number of layers to plot when ``y=None`` (default 16).
        main: Title string or list of strings (one per layer).
        ax: Existing matplotlib ``Axes`` to plot into (single-layer only).
        figsize: Figure size ``(width, height)`` in inches.
        **kwargs: Additional keyword arguments forwarded to
            ``matplotlib.pyplot.subplots``.

    Returns:
        * Single-layer: the ``matplotlib.axes.Axes`` used.
        * Multi-layer: 2-D ``numpy.ndarray`` of ``Axes``.
    """
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import builtins

    if not r.hasValues:
        warnings.warn("plot: SpatRaster has no cell values", stacklevel=2)

    nl_total = r.nlyr()
    lyr_names = list(r.names)

    # ── resolve layer selection (0-based) ────────────────────────────────────
    if y is None:
        n_plot = min(nl_total, maxnl)
        lyrs_0based = list(builtins.range(n_plot))
    elif isinstance(y, str):
        lyrs_0based = [lyr_names.index(y)]
    elif isinstance(y, (list, tuple)) and y and isinstance(y[0], str):
        lyrs_0based = [lyr_names.index(n) for n in y]
    elif isinstance(y, (list, tuple)):
        lyrs_0based = [int(v) for v in y]
    else:
        lyrs_0based = [int(y)]

    def _norm_lyr(i: int) -> int:
        j = int(i)
        if j < 0:
            j = nl_total + j
        if j < 0 or j >= nl_total:
            raise IndexError(f"layer index {i!r} out of range for nlyr={nl_total}")
        return j

    lyrs_0based = [_norm_lyr(i) for i in lyrs_0based]

    # ── colour palette ────────────────────────────────────────────────────────
    if col is None:
        palette = _default_palette(255)
    elif isinstance(col, str):
        cmap_obj = cm.get_cmap(col, 255)
        import matplotlib.colors as mc
        palette = [mc.to_hex(cmap_obj(i)) for i in builtins.range(255)]
    else:
        palette = list(col)

    # ── breaks imply interval type ────────────────────────────────────────────
    effective_type = type
    if breaks is not None and type is None:
        effective_type = "interval"

    # ── determine per-layer types ─────────────────────────────────────────────
    def _layer_type(lyr0: int) -> str:
        if effective_type is not None:
            return effective_type
        return _infer_type(r, lyr0)

    # ── single layer ──────────────────────────────────────────────────────────
    if len(lyrs_0based) == 1:
        lyr0 = lyrs_0based[0]
        ltype = _layer_type(lyr0)
        title_str = (main[0] if isinstance(main, (list, tuple)) else main) or lyr_names[lyr0]

        if ax is None:
            fig, ax_ = plt.subplots(1, 1, figsize=figsize or (6, 5), **kwargs)
        else:
            ax_ = ax

        _plot_one_layer(
            r, lyr0, ax_, palette, ltype,
            zlim=zlim, clamp=clamp, breaks=breaks,
            levels=levels, legend=legend, na_color=na_color,
            axes=axes, title=title_str,
            max_cell=maxcell, alpha=alpha, smooth=smooth,
        )
        return ax_

    # ── multi-layer ───────────────────────────────────────────────────────────
    n_plot = len(lyrs_0based)
    nrows, ncols = _get_nrnc(nr, nc, n_plot)

    if main is None:
        titles = [lyr_names[i] for i in lyrs_0based]
    elif isinstance(main, str):
        titles = [main] * n_plot
    else:
        titles = list(main) + [""] * max(0, n_plot - len(main))

    fig, axes_arr = plt.subplots(
        nrows, ncols,
        figsize=figsize or (4 * ncols, 3.5 * nrows),
        **kwargs,
    )

    axes_flat = np.array(axes_arr).flatten()

    for i, lyr0 in enumerate(lyrs_0based):
        ax_i = axes_flat[i]
        ltype = _layer_type(lyr0)
        _plot_one_layer(
            r, lyr0, ax_i, palette, ltype,
            zlim=zlim, clamp=clamp, breaks=breaks,
            levels=levels, legend=legend, na_color=na_color,
            axes=axes, title=titles[i],
            max_cell=maxcell // n_plot, alpha=alpha, smooth=smooth,
        )

    # Hide unused subplots
    for j in builtins.range(i + 1, len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.tight_layout()
    return np.array(axes_arr)


def plot_rgb(
    r: SpatRaster,
    red: int = 0,
    green: int = 1,
    blue: int = 2,
    alpha_band: Optional[int] = None,
    scale: float = 255.0,
    stretch: Optional[str] = None,
    smooth: bool = True,
    na_color: str = "white",
    axes: bool = False,
    ax: Optional[Any] = None,
    figsize: Optional[Tuple[float, float]] = None,
    **kwargs: Any,
) -> Any:
    """
    Plot a composite colour (RGB or RGBA) image from a multi-layer SpatRaster.

    Args:
        r: SpatRaster with at least three layers.
        red: 0-based layer index for the red channel.
        green: 0-based layer index for the green channel.
        blue: 0-based layer index for the blue channel.
        alpha_band: Optional 0-based layer index for the alpha (transparency)
            channel.
        scale: Maximum value for normalisation.  Use 255 for 8-bit data,
            65535 for 16-bit data, or 1 if bands are already in [0, 1].
        stretch: Contrast stretch method.  ``"lin"`` applies a linear
            2%–98% percentile stretch; ``None`` uses simple min-max
            normalisation to *scale*.
        smooth: If True (default), apply bilinear interpolation when
            rendering.
        na_color: Colour used for cells where any band is NA.
        axes: If True, draw coordinate axis ticks and labels.
        ax: Existing matplotlib ``Axes`` to plot into.
        figsize: Figure size ``(width, height)`` in inches.
        **kwargs: Additional keyword arguments forwarded to
            ``matplotlib.pyplot.subplots``.

    Returns:
        The ``matplotlib.axes.Axes`` used for the plot.
    """
    import matplotlib.pyplot as plt

    rgb_bands = [red, green, blue]
    if alpha_band is not None:
        rgb_bands.append(alpha_band)

    rgba = _rgb_image(r, rgb_bands, scale=scale, stretch=stretch, na_color=na_color)

    if ax is None:
        fig, ax_ = plt.subplots(1, 1, figsize=figsize or (6, 5), **kwargs)
    else:
        ax_ = ax

    ext_v = r.extent.vector
    lonlat = r.isLonLat()
    interp = "bilinear" if smooth else "nearest"

    _setup_axes(ax_, ext_v, axes, lonlat)
    ax_.imshow(rgba, extent=ext_v, origin="upper",
               interpolation=interp, aspect="auto")
    return ax_
