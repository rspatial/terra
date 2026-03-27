"""
write.py — write raster and vector data to files.
"""
from __future__ import annotations
import os
from typing import List, Optional, Union

from ._terra import SpatRaster, SpatVector, SpatOptions
from ._helpers import messages, spatoptions

_cpp_vect_write = SpatVector.write  # captured before monkey-patching


def _opt() -> SpatOptions:
    return SpatOptions()


# ---------------------------------------------------------------------------
# SpatRaster write
# ---------------------------------------------------------------------------

def write_raster(
    x: SpatRaster,
    filename: str,
    overwrite: bool = False,
    filetype: Optional[str] = None,
    datatype: str = "FLT4S",
    gdal: Optional[List[str]] = None,
    **kwargs,
) -> SpatRaster:
    """
    Write a SpatRaster to a file.

    Parameters
    ----------
    x : SpatRaster
    filename : str
        Output path.  The format is inferred from the extension
        (e.g. ``".tif"`` → GeoTIFF, ``".nc"`` → NetCDF).
    overwrite : bool
        If False (default), raise an error if *filename* already exists.
    filetype : str, optional
        GDAL driver name (e.g. ``"GTiff"``).  Auto-detected if None.
    datatype : str
        Output data type: ``"FLT4S"`` (float32, default), ``"FLT8S"``
        (float64), ``"INT2S"`` (int16), ``"INT4S"`` (int32), ``"INT1U"``
        (uint8), etc.
    gdal : list of str, optional
        GDAL creation options (e.g. ``["COMPRESS=LZW"]``).
    **kwargs
        Additional options forwarded to SpatOptions.

    Returns
    -------
    SpatRaster  (re-opened from the written file)
    """
    filename = os.path.expanduser(filename.strip())
    if not filename:
        raise ValueError("provide a filename")
    opt = spatoptions(filename, overwrite)
    opt.datatype = datatype
    if filetype:
        opt.filetype = filetype
    if gdal:
        opt.gdal_options = gdal
    xc = x.writeRaster(opt)
    messages(xc, "writeRaster")
    from .rast import rast
    return rast(filename)


def write_start(
    x: SpatRaster,
    filename: str,
    overwrite: bool = False,
    n: int = 4,
    sources: Optional[List[str]] = None,
) -> dict:
    """
    Open a file for block-wise writing.

    Parameters
    ----------
    x : SpatRaster
    filename : str
    overwrite : bool
    n : int
        Number of block copies to buffer in memory.
    sources : list of str, optional

    Returns
    -------
    dict with keys ``"n"``, ``"row"``, ``"nrows"``.
    """
    filename = os.path.expanduser(filename.strip())
    opt = spatoptions(filename, overwrite)
    opt.ncopies = n
    if sources is None:
        sources = []
    ok = x.writeStart(opt, list(set(sources)))
    messages(x, "writeStart")
    b = x.getBlockSizeWrite()
    b["row"] = [r + 1 for r in b["row"]]
    return b


def write_values(
    x: SpatRaster,
    v: List[float],
    start: int,
    nrows: int,
) -> bool:
    """
    Write a block of values to a file opened with write_start().

    Parameters
    ----------
    x : SpatRaster
    v : list of float
    start : int
        Starting row (1-based).
    nrows : int
        Number of rows in this block.

    Returns
    -------
    bool
    """
    ok = x.writeValues(list(v), start - 1, nrows)
    messages(x, "writeValues")
    return bool(ok)


def write_stop(x: SpatRaster) -> SpatRaster:
    """
    Finalise a block-wise write and close the file.

    Parameters
    ----------
    x : SpatRaster

    Returns
    -------
    SpatRaster  (re-opened from the written file)
    """
    x.writeStop()
    messages(x, "writeStop")
    src = list(x.filenames())
    if src and src[0]:
        from .rast import rast
        return rast(src[0])
    return x


def blocks(x: SpatRaster, n: int = 4) -> dict:
    """
    Return the block structure for writing *x*.

    Parameters
    ----------
    x : SpatRaster
    n : int
        Number of copies to buffer.

    Returns
    -------
    dict with keys ``"n"``, ``"row"`` (1-based), ``"nrows"``.
    """
    opt = spatoptions()
    opt.ncopies = n
    b = x.getBlockSizeR(opt)
    b["row"] = [r + 1 for r in b["row"]]
    return b


# ---------------------------------------------------------------------------
# SpatVector write
# ---------------------------------------------------------------------------

_EXT_TO_FILETYPE = {
    "shp": "ESRI Shapefile",
    "shz": "ESRI Shapefile",
    "gpkg": "GPKG",
    "gdb": "OpenFileGDB",
    "gml": "GML",
    "json": "GeoJSON",
    "geojson": "GeoJSON",
    "cdf": "netCDF",
    "svg": "SVG",
    "kml": "KML",
    "vct": "Idrisi",
    "tab": "MapInfo File",
}


def _guess_filetype(filename: str) -> str:
    ext = os.path.splitext(filename)[1].lstrip(".").lower()
    ft = _EXT_TO_FILETYPE.get(ext)
    if ft is None:
        raise ValueError(f"cannot guess filetype from filename {filename!r}")
    return ft


def write_vector(
    x: SpatVector,
    filename: str,
    filetype: Optional[str] = None,
    layer: Optional[str] = None,
    insert: bool = False,
    overwrite: bool = False,
    options: str = "ENCODING=UTF-8",
) -> bool:
    """
    Write a SpatVector to a file.

    Parameters
    ----------
    x : SpatVector
    filename : str
        Output path.
    filetype : str, optional
        OGR driver name.  Auto-detected from extension if None.
    layer : str, optional
        Layer name inside the file.  Defaults to the base filename.
    insert : bool
        Append to an existing file instead of replacing.
    overwrite : bool
        Overwrite an existing layer.
    options : str
        OGR creation options string.

    Returns
    -------
    bool
    """
    filename = os.path.expanduser(filename.strip())
    if not filename:
        raise ValueError("provide a filename")
    if filetype is None:
        filetype = _guess_filetype(filename)
    if layer is None:
        layer = os.path.splitext(os.path.basename(filename))[0].strip()

    # Truncate field names for Shapefiles (max 10 chars)
    if filetype == "ESRI Shapefile":
        nms = list(x.names)
        truncated = [n[:10] for n in nms]
        if truncated != nms:
            xc = x.deepcopy()
            xc.names = truncated
            x = xc

    ok = _cpp_vect_write(x, filename, layer, filetype, insert, overwrite, options)
    messages(x, "writeVector")
    return bool(ok)
