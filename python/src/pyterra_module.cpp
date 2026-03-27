// pyterra_module.cpp  — pybind11 bindings for the terra C++ core
// Parallel to src/RcppModule.cpp but targeting Python instead of R.
//
// Build with -Dstandalone so spatBase.h does NOT define useRcpp.
// useGDAL remains active (no -Dnogdal).

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       // std::vector<T> <-> Python list auto-conversion
#include <cstdint>

#include "spatRasterMultiple.h"
#include "spatGraph.h"
#include <memory>
#include "NA.h"
#include "spatTime.h"

// GDAL headers for module-level initialisation
#include "gdal_priv.h"
#include "ogr_api.h"
#include "cpl_conv.h"

namespace py = pybind11;


// ── Bridge helpers ──────────────────────────────────────────────────────────
// These replace the Rcpp-returning functions in RcppModule.cpp.

// SpatDataFrame → Python dict  {colname: list}
static py::dict getDataFrame(SpatDataFrame* v) {
    py::dict out;
    size_t n = v->ncol();
    if (n == 0) return out;

    auto longNA   = NA<long>::value;
    std::string stringNA = v->NAS;
    SpatTime_t timeNA = NA<SpatTime_t>::value;

    std::vector<std::string> nms   = v->names;
    std::vector<size_t>      itype = v->itype;

    for (size_t i = 0; i < n; i++) {
        const std::string& nm = nms[i];
        if (itype[i] == 0) {
            out[py::str(nm)] = v->getD(i);                     // double
        } else if (itype[i] == 1) {
            std::vector<long> ints = v->getI(i);
            py::list col;
            for (auto x : ints)
                col.append(x == longNA ? py::object(py::none()) : py::object(py::int_(x)));
            out[py::str(nm)] = col;
        } else if (itype[i] == 2) {
            std::vector<std::string> s = v->getS(i);
            py::list col;
            for (auto& x : s)
                col.append(x == stringNA ? py::object(py::none()) : py::object(py::str(x)));
            out[py::str(nm)] = col;
        } else if (itype[i] == 3) {
            std::vector<int8_t> b = v->getB(i);
            py::list col;
            for (auto x : b)
                col.append(x > 1 ? py::object(py::none()) : py::object(py::bool_(x != 0)));
            out[py::str(nm)] = col;
        } else if (itype[i] == 4) {
            SpatTime_v tx = v->getT(i);
            py::list col;
            for (auto x : tx.x)
                col.append(x == timeNA ? py::object(py::none()) : py::object(py::float_(static_cast<double>(x))));
            out[py::str(nm)] = col;
        } else if (itype[i] == 5) {
            out[py::str(nm)] = py::cast(v->getF(i));           // SpatFactor (registered below)
        }
    }
    return out;
}

static py::dict getVectorAttributes(SpatVector* v) {
    SpatDataFrame df = v->df;
    return getDataFrame(&df);
}

// geometry as a flat dict-of-lists  {id, part, x, y, hole}
static py::dict get_geometryDF(SpatVector* v) {
    SpatDataFrame df = v->getGeometryDF();
    py::dict out;
    out["id"]   = df.iv[0];
    out["part"] = df.iv[1];
    out["x"]    = df.dv[0];
    out["y"]    = df.dv[1];
    out["hole"] = df.iv[2];
    return out;
}

// geometry as nested list-of-list-of-dicts: [geom][part/hole]{x,y}
static py::list get_geometryList(SpatVector* v,
                                  const std::string xnm,
                                  const std::string ynm) {
    py::list out;
    size_t ni = v->nrow();
    for (size_t i = 0; i < ni; i++) {
        SpatGeom g = v->getGeom(i);
        size_t nj = g.size();
        py::list geom;
        for (size_t j = 0; j < nj; j++) {
            SpatPart p = g.getPart(j);
            py::list parts;
            py::dict m;
            m[py::str(xnm)] = p.x;
            m[py::str(ynm)] = p.y;
            parts.append(m);
            for (size_t k = 0; k < p.nHoles(); k++) {
                SpatHole h = p.getHole(k);
                py::dict mh;
                mh[py::str(xnm)] = h.x;
                mh[py::str(ynm)] = h.y;
                parts.append(mh);
            }
            geom.append(parts);
        }
        out.append(geom);
    }
    return out;
}

static py::dict getBlockSizeR(SpatRaster* r, SpatOptions* opt) {
    BlockSize bs = r->getBlockSize(*opt);
    py::dict d;
    d["row"]   = bs.row;
    d["nrows"] = bs.nrows;
    d["n"]     = bs.n;
    return d;
}

static py::dict getBlockSizeWrite(SpatRaster* r) {
    BlockSize bs = r->bs;
    py::dict d;
    d["row"]   = bs.row;
    d["nrows"] = bs.nrows;
    d["n"]     = bs.n;
    return d;
}

// Accept a Python list of bytes objects (WKB geometries)
static bool addWKB(SpatVector* x, py::list wkb) {
    size_t n = py::len(wkb);
    // Keep py::bytes objects alive while addRawGeoms runs
    std::vector<py::bytes> keep(n);
    std::vector<unsigned char*> raw(n);
    std::vector<size_t> sizes(n);
    for (size_t i = 0; i < n; i++) {
        keep[i] = wkb[i].cast<py::bytes>();
        char* buf; Py_ssize_t len;
        PyBytes_AsStringAndSize(keep[i].ptr(), &buf, &len);
        raw[i]   = reinterpret_cast<unsigned char*>(buf);
        sizes[i] = static_cast<size_t>(len);
    }
    return x->addRawGeoms(raw, sizes);
}


// ── Module definition ───────────────────────────────────────────────────────
PYBIND11_MODULE(_pyterra, m) {
    m.doc() = "pyterra: Python bindings for the terra spatial C++ library";

    // Initialise GDAL/OGR once when the module is imported
    GDALAllRegister();
    OGRRegisterAll();
    CPLSetConfigOption("GDAL_MAX_BAND_COUNT", "9999999");
    CPLSetConfigOption("OGR_CT_FORCE_TRADITIONAL_GIS_ORDER", "YES");


    // ── SpatTime_v ──────────────────────────────────────────────────────────
    py::class_<SpatTime_v>(m, "SpatTime_v")
        .def(py::init<>())
        .def_readwrite("step", &SpatTime_v::step)
        .def_readwrite("zone", &SpatTime_v::zone)
        .def_readwrite("x",    &SpatTime_v::x)
    ;


    // ── SpatFactor ──────────────────────────────────────────────────────────
    py::class_<SpatFactor>(m, "SpatFactor")
        .def(py::init<>())
        .def(py::init<std::vector<size_t>, std::vector<std::string>, bool>())
        .def_readwrite("values",  &SpatFactor::v)
        .def_readwrite("labels",  &SpatFactor::labels)
        .def_readwrite("ordered", &SpatFactor::ordered)
    ;


    // ── SpatSRS ─────────────────────────────────────────────────────────────
    py::class_<SpatSRS>(m, "SpatSRS")
        .def(py::init<>())
        .def("set",       &SpatSRS::set)
        .def("is_lonlat", &SpatSRS::is_lonlat)
        .def("to_meter",  &SpatSRS::to_meter)
    ;


    // ── SpatExtent ──────────────────────────────────────────────────────────
    py::class_<SpatExtent>(m, "SpatExtent")
        .def(py::init<>())
        .def(py::init<double, double, double, double>())
        .def("deepcopy",       &SpatExtent::deepCopy)
        .def_property_readonly("vector",       &SpatExtent::asVector)
        .def_property_readonly("valid",        &SpatExtent::valid)
        .def_property_readonly("valid_notempty",&SpatExtent::valid_notempty)
        .def_property_readonly("empty",        &SpatExtent::empty)
        .def("align",          &SpatExtent::align)
        .def("intersect",      &SpatExtent::intersect)
        .def("as_points",      &SpatExtent::asPoints)
        .def("ceil",           &SpatExtent::ceil)
        .def("compare",        &SpatExtent::compare)
        .def("floor",          &SpatExtent::floor)
        .def("round",          &SpatExtent::round)
        .def("union",          &SpatExtent::unite)
        .def("sampleRandom",   &SpatExtent::sampleRandom)
        .def("sampleRegular",  &SpatExtent::sampleRegular)
        .def("sample",         &SpatExtent::test_sample)
    ;


    // ── SpatMessages ────────────────────────────────────────────────────────
    py::class_<SpatMessages>(m, "SpatMessages")
        .def(py::init<>())
        .def_readwrite("has_error",   &SpatMessages::has_error)
        .def_readwrite("has_warning", &SpatMessages::has_warning)
    ;


    // ── SpatOptions ─────────────────────────────────────────────────────────
    py::class_<SpatOptions>(m, "SpatOptions")
        .def(py::init<>())
        .def("deepcopy", &SpatOptions::deepCopy)
        .def_readwrite("parallel",      &SpatOptions::parallel)
        .def_readwrite("metadata",      &SpatOptions::tags)
        .def_property("tempdir",   &SpatOptions::get_tempdir, &SpatOptions::set_tempdir)
        .def_property("memfrac",   &SpatOptions::get_memfrac, &SpatOptions::set_memfrac)
        .def_property("memmax",    &SpatOptions::get_memmax,  &SpatOptions::set_memmax)
        .def_property("memmin",    &SpatOptions::get_memmin,  &SpatOptions::set_memmin)
        .def_property("tolerance", &SpatOptions::get_tolerance, &SpatOptions::set_tolerance)
        .def_property("filenames", &SpatOptions::get_filenames, &SpatOptions::set_filenames)
        .def_property("filetype",  &SpatOptions::get_filetype,  &SpatOptions::set_filetype)
        .def_property("datatype",  &SpatOptions::get_datatype,  &SpatOptions::set_datatype)
        .def_property("verbose",   &SpatOptions::get_verbose,   &SpatOptions::set_verbose)
        .def_property("NAflag",    &SpatOptions::get_NAflag,    &SpatOptions::set_NAflag)
        .def_property("statistics",&SpatOptions::get_statistics,&SpatOptions::set_statistics)
        .def_property("overwrite", &SpatOptions::get_overwrite, &SpatOptions::set_overwrite)
        .def_readwrite("datatype_set", &SpatOptions::datatype_set)
        .def_readwrite("threads",      &SpatOptions::threads)
        .def_property("progress",  &SpatOptions::get_progress,  &SpatOptions::set_progress)
        .def_readwrite("progressbar",  &SpatOptions::progressbar)
        .def_property("ncopies",   &SpatOptions::get_ncopies,   &SpatOptions::set_ncopies)
        .def_property("def_filetype",&SpatOptions::get_def_filetype,&SpatOptions::set_def_filetype)
        .def_property("def_datatype",&SpatOptions::get_def_datatype,&SpatOptions::set_def_datatype)
        .def_property("def_verbose", &SpatOptions::get_def_verbose, &SpatOptions::set_def_verbose)
        .def_property("todisk",    &SpatOptions::get_todisk,    &SpatOptions::set_todisk)
        .def_readwrite("messages",     &SpatOptions::msg)
        .def_readwrite("gdal_options", &SpatOptions::gdal_options)
        .def_readwrite("tmpfile",      &SpatOptions::tmpfile)
        .def_readwrite("names",        &SpatOptions::names)
        .def_property("steps",     &SpatOptions::get_steps,     &SpatOptions::set_steps)
        .def_property("scale",     &SpatOptions::get_scale,     &SpatOptions::set_scale)
        .def_property("offset",    &SpatOptions::get_offset,    &SpatOptions::set_offset)
        .def("has_error",   &SpatOptions::hasError)
        .def("has_warning", &SpatOptions::hasWarning)
        .def("getError",    &SpatOptions::getError)
        .def("getWarnings", &SpatOptions::getWarnings)
    ;


    // ── SpatDataFrame ───────────────────────────────────────────────────────
    py::class_<SpatDataFrame>(m, "SpatDataFrame")
        .def(py::init<>())
        .def_readonly("itype",  &SpatDataFrame::itype)
        .def_readonly("iplace", &SpatDataFrame::iplace)
        .def_property("names",
            &SpatDataFrame::get_names,
            &SpatDataFrame::set_names)
        .def_property("nrow",
            (size_t (SpatDataFrame::*)() const)   &SpatDataFrame::nrow,
            (void   (SpatDataFrame::*)(size_t))    &SpatDataFrame::resize_rows)
        .def_property("ncol",
            (size_t (SpatDataFrame::*)() const)   &SpatDataFrame::ncol,
            (void   (SpatDataFrame::*)(size_t))    &SpatDataFrame::resize_cols)
        .def("has_error",   &SpatDataFrame::hasError)
        .def("has_warning", &SpatDataFrame::hasWarning)
        .def("getWarnings", &SpatDataFrame::getWarnings)
        .def("getError",    &SpatDataFrame::getError)
        .def("add_column_double",
            (bool (SpatDataFrame::*)(std::vector<double>, std::string))(&SpatDataFrame::add_column))
        .def("add_column_long",
            (bool (SpatDataFrame::*)(std::vector<long>, std::string))(&SpatDataFrame::add_column))
        .def("add_column_string",
            (bool (SpatDataFrame::*)(std::vector<std::string>, std::string))(&SpatDataFrame::add_column))
        .def("add_column_factor",
            (bool (SpatDataFrame::*)(SpatFactor, std::string))(&SpatDataFrame::add_column))
        .def("add_column_bool",
            (bool (SpatDataFrame::*)(std::vector<int>, std::string))(&SpatDataFrame::add_column_bool))
        .def("add_column_time",   &SpatDataFrame::add_column_time)
        .def("remove_column",
            (bool (SpatDataFrame::*)(std::string))(&SpatDataFrame::remove_column))
        .def("remove_column_i",
            (bool (SpatDataFrame::*)(int))(&SpatDataFrame::remove_column))
        .def("get_datatypes",  &SpatDataFrame::get_datatypes)
        .def("get_timezones",  &SpatDataFrame::get_timezones)
        .def("get_timesteps",  &SpatDataFrame::get_timesteps)
        .def("subset_rows",
            (SpatDataFrame (SpatDataFrame::*)(std::vector<size_t>))(&SpatDataFrame::subset_rows))
        .def("subset_cols",
            (SpatDataFrame (SpatDataFrame::*)(std::vector<size_t>))(&SpatDataFrame::subset_cols))
        .def("remove_rows",  &SpatDataFrame::remove_rows)
        .def("cbind",        &SpatDataFrame::cbind)
        .def("rbind",        &SpatDataFrame::rbind)
        .def("values",       &getDataFrame)      // returns Python dict
        .def("unique",       &SpatDataFrame::unique)
        .def("write",        &SpatDataFrame::write_dbf)
        .def_readwrite("messages", &SpatDataFrame::msg)
        .def("strwidth",     &SpatDataFrame::strwidth)
    ;


    // ── SpatCategories ──────────────────────────────────────────────────────
    py::class_<SpatCategories>(m, "SpatCategories")
        .def(py::init<>())
        .def_readonly("df",     &SpatCategories::d)
        .def_readwrite("index", &SpatCategories::index)
        .def("combine",         &SpatCategories::combine)
    ;


    // ── SpatVectorCollection ────────────────────────────────────────────────
    py::class_<SpatVectorCollection>(m, "SpatVectorCollection")
        .def(py::init<>())
        .def(py::init<std::string, std::string, std::string, std::string,
                      std::vector<double>, SpatVector>())
        .def("deepcopy",    &SpatVectorCollection::deepCopy)
        .def("size",        &SpatVectorCollection::size)
        .def("get",         &SpatVectorCollection::get)
        .def("push_back",   &SpatVectorCollection::push_back)
        .def("subset",      &SpatVectorCollection::subset)
        .def("replace",     &SpatVectorCollection::replace)
        .def("append",      &SpatVectorCollection::append)
        .def("has_error",   &SpatVectorCollection::hasError)
        .def("has_warning", &SpatVectorCollection::hasWarning)
        .def("getWarnings", &SpatVectorCollection::getWarnings)
        .def("getError",    &SpatVectorCollection::getError)
        .def("from_hex_col",&SpatVectorCollection::from_hex_col)
        .def("setNames",    &SpatVectorCollection::setNames)
        .def_property_readonly("names", &SpatVectorCollection::getNames)
    ;


    // ── SpatVector ──────────────────────────────────────────────────────────
    py::class_<SpatVector>(m, "SpatVector")
        .def(py::init<>())
        .def(py::init<SpatExtent, std::string>())
        .def(py::init<std::vector<std::string>>())

        .def("addWKB",       &addWKB)
        .def("deepcopy",     &SpatVector::deepCopy)
        .def("wkt",          &SpatVector::wkt)
        .def("wkb",          &SpatVector::wkb)
        .def("wkb_raw",      &SpatVector::wkb_raw)
        .def("hex",          &SpatVector::hex)
        .def("from_hex",     &SpatVector::from_hex)
        .def("make_nodes",   &SpatVector::make_nodes)
        .def("boundary",     &SpatVector::boundary)
        .def("polygonize",   &SpatVector::polygonize)
        .def("normalize",    &SpatVector::normalize)
        .def("normalize_longitude", &SpatVector::normalize_longitude)
        .def("rotate_longitude",    &SpatVector::rotate_longitude)
        .def("line_merge",   &SpatVector::line_merge)
        .def("simplify",     &SpatVector::simplify)
        .def("thin",         &SpatVector::thin)
        .def("shared_paths",
            (SpatVector (SpatVector::*)(bool))(&SpatVector::shared_paths))
        .def("shared_paths2",
            (SpatVector (SpatVector::*)(SpatVector, bool))(&SpatVector::shared_paths))
        .def("snap",         &SpatVector::snap)
        .def("snapto",       &SpatVector::snapto)
        .def("spatial_index_2d",     &SpatVector::index_2d)
        .def("spatial_index_sparse", &SpatVector::index_sparse)

        .def_readonly("is_proxy",    &SpatVector::is_proxy)
        .def_readonly("read_query",  &SpatVector::read_query)
        .def_readonly("read_extent", &SpatVector::read_extent)
        .def_readonly("geom_count",  &SpatVector::geom_count)
        .def_readonly("source",      &SpatVector::source)
        .def_readonly("layer",       &SpatVector::source_layer)
        .def_readonly("df",          &SpatVector::df)

        .def("has_error",   &SpatVector::hasError)
        .def("has_warning", &SpatVector::hasWarning)
        .def("getError",    &SpatVector::getError)
        .def("getWarnings", &SpatVector::getWarnings)

        .def("coordinates",     &SpatVector::coordinates)
        .def("get_geometry",    &SpatVector::getGeometry)
        .def("get_geometryDF",  &get_geometryDF)
        .def("get_geometryList",&get_geometryList)
        .def("linesList",       &SpatVector::linesList)
        .def("polygonsList",    &SpatVector::polygonsList)
        .def("linesNA",         &SpatVector::linesNA)

        .def("add_column_empty",
            (void (SpatVector::*)(unsigned, std::string))(&SpatVector::add_column))
        .def("add_column_double",
            (bool (SpatVector::*)(std::vector<double>, std::string))(&SpatVector::add_column))
        .def("add_column_long",
            (bool (SpatVector::*)(std::vector<long>, std::string))(&SpatVector::add_column))
        .def("add_column_string",
            (bool (SpatVector::*)(std::vector<std::string>, std::string))(&SpatVector::add_column))
        .def("add_column_factor", &SpatVector::add_column_factor)
        .def("add_column_bool",   &SpatVector::add_column_bool)
        .def("add_column_time",   &SpatVector::add_column_time)
        .def("remove_column",
            (bool (SpatVector::*)(std::string))(&SpatVector::remove_column))
        .def("remove_column_i",
            (bool (SpatVector::*)(int))(&SpatVector::remove_column))
        .def("remove_df",   &SpatVector::remove_df)
        .def("set_df",      &SpatVector::set_df)
        .def("get_datatypes",&SpatVector::get_datatypes)

        .def("set_holes",    &SpatVector::set_holes)
        .def("get_holes",    &SpatVector::get_holes)
        .def("remove_holes", &SpatVector::remove_holes)
        .def("append",       &SpatVector::append)
        .def("cbind",        &SpatVector::cbind)

        .def("area",         &SpatVector::area)
        .def("as_lines",     &SpatVector::as_lines)
        .def("as_points",    &SpatVector::as_points)
        .def("couldBeLonLat",&SpatVector::could_be_lonlat)
        .def("get_crs",      &SpatVector::getSRS)
        .def("set_crs",
            (bool (SpatVector::*)(std::string))(&SpatVector::setSRS))

        .def("distance_self",
            (std::vector<double> (SpatVector::*)(bool, std::string, const std::string, bool, SpatOptions&))(&SpatVector::distance))
        .def("distance_other",
            (std::vector<double> (SpatVector::*)(SpatVector, bool, std::string, const std::string, bool, SpatOptions&))(&SpatVector::distance))
        .def("point_distance", &SpatVector::pointdistance)

        .def("extent",       &SpatVector::getExtent)
        .def("getDF",        &getVectorAttributes)    // returns Python dict
        .def("getGeometryWKT",&SpatVector::getGeometryWKT)
        .def("isLonLat",     &SpatVector::is_lonlat)
        .def("length",       &SpatVector::length)
        .def("nsegments",    &SpatVector::nseg)
        .def_readwrite("messages", &SpatVector::msg)
        .def_property("names",
            &SpatVector::get_names,
            &SpatVector::set_names)
        .def("nrow",         &SpatVector::nrow)
        .def("ncol",         &SpatVector::ncol)
        .def("project",      &SpatVector::project)
        .def("project_xy",   &SpatVector::project_xy)
        .def("read",         &SpatVector::read)
        .def("setGeometry",  &SpatVector::setGeometry)
        .def("setPointsXY",  &SpatVector::setPointsGeometry)
        .def("setPointsDF",  &SpatVector::setPointsDF)
        .def("setLinesStartEnd",&SpatVector::setLinesStartEnd)
        .def("size",         &SpatVector::size)
        .def("subset_cols",
            (SpatVector (SpatVector::*)(std::vector<long>))(&SpatVector::subset_cols))
        .def("subset_rows",
            (SpatVector (SpatVector::*)(std::vector<long>))(&SpatVector::subset_rows))
        .def("remove_rows",  &SpatVector::remove_rows)
        .def("type",         &SpatVector::type)
        .def("multipoint",   &SpatVector::is_multipoint)
        .def("naGeoms",      &SpatVector::naGeoms)
        .def("nullGeoms",    &SpatVector::nullGeoms)
        .def("write",        &SpatVector::write)
        .def("delete_layers",&SpatVector::delete_layers)
        .def("layer_names",  &SpatVector::layer_names)

        .def("geos_isvalid",     &SpatVector::geos_isvalid)
        .def("geos_isvalid_msg", &SpatVector::geos_isvalid_msg)

        .def("aggregate",
            (SpatVector (SpatVector::*)(std::string, bool))(&SpatVector::aggregate))
        .def("aggregate_nofield",
            (SpatVector (SpatVector::*)(bool))(&SpatVector::aggregate))
        .def("disaggregate", &SpatVector::disaggregate)
        .def("buffer",       &SpatVector::buffer)
        .def("centroid",     &SpatVector::centroid)
        .def("point_on_surface", &SpatVector::point_on_surface)
        .def("make_valid2",  &SpatVector::make_valid2)
        .def("flip",         &SpatVector::flip)
        .def("transpose",    &SpatVector::transpose)
        .def("shift",        &SpatVector::shift)
        .def("rescale",      &SpatVector::rescale)
        .def("rotate",       &SpatVector::rotate)
        .def("erase_agg",    &SpatVector::erase_agg)
        .def("erase",
            (SpatVector (SpatVector::*)(SpatVector))(&SpatVector::erase))
        .def("erase_self",
            (SpatVector (SpatVector::*)(bool))(&SpatVector::erase))
        .def("elongate",     &SpatVector::elongate)
        .def("gaps",         &SpatVector::gaps)
        .def("symdif",       &SpatVector::symdif)
        .def("cover",        &SpatVector::cover)
        .def("union",
            (SpatVector (SpatVector::*)(SpatVector))(&SpatVector::unite))
        .def("union_self",
            (SpatVector (SpatVector::*)())(&SpatVector::unite))
        .def("union_unary",  &SpatVector::unaryunion)
        .def("intersect",    &SpatVector::intersect)
        .def("delaunay",     &SpatVector::delaunay)
        .def("voronoi",      &SpatVector::voronoi)
        .def("hull",         &SpatVector::hull)
        .def("width",        &SpatVector::width)
        .def("clearance",    &SpatVector::clearance)
        .def("mask",         &SpatVector::mask)
        .def("is_related",   &SpatVector::is_related)
        .def("related_between",
            (std::vector<std::vector<double>> (SpatVector::*)(SpatVector, std::string, bool))(&SpatVector::which_relate))
        .def("related_within",
            (std::vector<std::vector<double>> (SpatVector::*)(std::string, bool))(&SpatVector::which_relate))
        .def("relate_within",
            (std::vector<int> (SpatVector::*)(std::string, bool))(&SpatVector::relate))
        .def("equals_between",
            (std::vector<size_t> (SpatVector::*)(SpatVector, double))(&SpatVector::equals_exact))
        .def("equals_within",
            (std::vector<size_t> (SpatVector::*)(bool, double))(&SpatVector::equals_exact))
        .def("crop_ext",
            (SpatVector (SpatVector::*)(SpatExtent, bool))(&SpatVector::crop))
        .def("crop_vct",
            (SpatVector (SpatVector::*)(SpatVector))(&SpatVector::crop))
        .def("near_between",
            (SpatVector (SpatVector::*)(SpatVector, bool, const std::string))(&SpatVector::nearest_point))
        .def("near_within",
            (SpatVector (SpatVector::*)(const std::string))(&SpatVector::nearest_point))
        .def("near_geom",    &SpatVector::nearest_geometry)
        .def("split",        &SpatVector::split)
        .def("sample",       &SpatVector::sample)
        .def("sampleGeom",   &SpatVector::sample_geom)
        .def("remove_duplicate_nodes", &SpatVector::remove_duplicate_nodes)
        .def("cross_dateline",         &SpatVector::cross_dateline)
        .def("fix_lonlat_overflow",    &SpatVector::fix_lonlat_overflow)
        .def("densify",      &SpatVector::densify)
        .def("round",        &SpatVector::round)
        .def("make_CCW",     &SpatVector::make_CCW)
    ;


    // ── SpatVectorProxy ─────────────────────────────────────────────────────
    py::class_<SpatVectorProxy>(m, "SpatVectorProxy")
        .def(py::init<>())
        .def_readwrite("v",  &SpatVectorProxy::v)
        .def("deepcopy",     &SpatVectorProxy::deepCopy)
    ;


    // ── SpatRaster ──────────────────────────────────────────────────────────
    py::class_<SpatRaster>(m, "SpatRaster")
        .def(py::init<>())
        .def(py::init<std::vector<std::string>, std::vector<int>,
                      std::vector<std::string>, bool,
                      std::vector<std::string>, std::vector<std::string>,
                      std::vector<int>, bool, bool,
                      std::vector<std::string>>())
        .def(py::init<std::vector<size_t>, std::vector<double>, std::string>())

        .def("has_error",   &SpatRaster::hasError)
        .def("has_warning", &SpatRaster::hasWarning)
        .def("getError",    &SpatRaster::getError)
        .def("getWarnings", &SpatRaster::getWarnings)
        .def("getMessage",  &SpatRaster::getMessage)

        .def("addTag",      &SpatRaster::addTag)
        .def("getTags",     &SpatRaster::getTags)
        .def("addLyrTags",  &SpatRaster::addLyrTags)
        .def("getLyrTags",  &SpatRaster::getLyrTags)
        .def("getAllFiles",  &SpatRaster::getAllFiles)
        .def("getFileBlocksize", &SpatRaster::getFileBlocksize)

        .def("sources_to_disk", &SpatRaster::sources_to_disk)
        .def("mem_needs",       &SpatRaster::mem_needs)
        .def("spatinit",        &SpatRaster::gdalogrproj_init)
        .def("addSource",       &SpatRaster::addSource)
        .def("replace",         &SpatRaster::replace)
        .def("combineSources",  &SpatRaster::combineSources)
        .def("compare_geom",    &SpatRaster::compare_geom)
        .def("couldBeLonLat",   &SpatRaster::could_be_lonlat)
        .def("deepcopy",        &SpatRaster::deepCopy)
        .def("hardcopy",        &SpatRaster::hardCopy)
        .def("get_crs",         &SpatRaster::getSRS)
        .def("set_crs",
            (bool (SpatRaster::*)(std::string))(&SpatRaster::setSRS))
        .def_property("extent",
            &SpatRaster::getExtent,
            (void (SpatRaster::*)(SpatExtent))&SpatRaster::setExtent)
        .def("is_rotated",  &SpatRaster::is_rotated)
        .def("is_flipped",  &SpatRaster::is_flipped)
        .def("setWindow",   &SpatRaster::setWindow)
        .def("removeWindow",&SpatRaster::removeWindow)
        .def("hasWindow",   &SpatRaster::hasWindow)
        .def("filenames",   &SpatRaster::filenames)

        .def_readwrite("rgb",    &SpatRaster::rgb)
        .def_readwrite("rgbtype",&SpatRaster::rgbtype)
        .def("setRGB",    &SpatRaster::setRGB)
        .def("removeRGB", &SpatRaster::removeRGB)
        .def("getRGB",    &SpatRaster::getRGB)

        .def("makeCategorical",  &SpatRaster::makeCategorical)
        .def("createCategories", &SpatRaster::createCategories)
        .def("hasCategories",    &SpatRaster::hasCategories)
        .def("getCategories",    &SpatRaster::getCategories)
        .def("setCategories",    &SpatRaster::setCategories)
        .def("removeCategories", &SpatRaster::removeCategories)
        .def("getLabels",        &SpatRaster::getLabels)
        .def("setLabels",        &SpatRaster::setLabels)
        .def("getCatIndex",      &SpatRaster::getCatIndex)
        .def("setCatIndex",      &SpatRaster::setCatIndex)
        .def("getScaleOffset",   &SpatRaster::getScaleOffset)
        .def("setScaleOffset",   &SpatRaster::setScaleOffset)
        .def("hasColors",        &SpatRaster::hasColors)
        .def("getColors",        &SpatRaster::getColors)
        .def("setColors",        &SpatRaster::setColors)
        .def("removeColors",     &SpatRaster::removeColors)

        .def_property_readonly("dataType",  &SpatRaster::dataType)
        .def("getDataType",  &SpatRaster::getDataType)
        .def("valueType",    &SpatRaster::getValueType)
        .def("setValueType", &SpatRaster::setValueType)
        .def_property_readonly("hasRange",  &SpatRaster::hasRange)
        .def_property_readonly("hasValues", &SpatRaster::hasValues)
        .def_property_readonly("inMemory",  &SpatRaster::inMemory)
        .def("isLonLat",        &SpatRaster::is_lonlat)
        .def("isGlobalLonLat",  &SpatRaster::is_global_lonlat)
        .def_property_readonly("names",     &SpatRaster::getNames)
        .def("get_sourcenames_long", &SpatRaster::getLongSourceNames)
        .def("set_sourcenames_long", &SpatRaster::setLongSourceNames)
        .def("get_sourcenames",      &SpatRaster::getSourceNames)
        .def("set_sourcenames",      &SpatRaster::setSourceNames)
        .def("setNAflag", &SpatRaster::setNAflag)
        .def("getNAflag", &SpatRaster::getNAflag)
        .def_property_readonly("isMD",     &SpatRaster::isMD)
        .def_property_readonly("hasUnit",  &SpatRaster::hasUnit)
        .def_property_readonly("hasDepth", &SpatRaster::hasDepth)
        .def_property_readonly("hasTime",  &SpatRaster::hasTime)
        .def_property_readonly("time",     &SpatRaster::getTime)
        .def_property_readonly("timestep", &SpatRaster::getTimeStep)
        .def_property_readonly("timezone", &SpatRaster::getTimeZone)

        .def("setTime",     &SpatRaster::setTime)
        .def("setDepth",    &SpatRaster::setDepth)
        .def("setUnit",     &SpatRaster::setUnit)
        .def("getUnit",     &SpatRaster::getUnit)
        .def("getDepth",    &SpatRaster::getDepth)

        .def("ncol",   &SpatRaster::ncol)
        .def("nrow",   &SpatRaster::nrow)
        .def("nlyr",   &SpatRaster::nlyr)
        .def("size",   &SpatRaster::size)
        .def("ncell",  &SpatRaster::ncell)
        .def("res",    &SpatRaster::resolution)
        .def("set_resolution", &SpatRaster::setResolution)

        .def("getBlockSize",  &getBlockSizeR)
        .def("getBlockSizeWrite", &getBlockSizeWrite)

        .def("readStart", &SpatRaster::readStart)
        .def("readStop",  &SpatRaster::readStop)
        .def("readAll",   &SpatRaster::readAll)
        .def("readValues", &SpatRaster::readValuesR)
        .def("getValues", &SpatRaster::getValues)

        .def("setValues",
            (bool (SpatRaster::*)(std::vector<double>&, SpatOptions&))(&SpatRaster::setValues))

        .def("writeStart", &SpatRaster::writeStart)
        .def("writeStop",  &SpatRaster::writeStop)
        .def("writeValues",&SpatRaster::writeValues)
        .def("writeValuesRect",&SpatRaster::writeValuesRect)

        .def("setRange",  &SpatRaster::setRange)
        .def_property_readonly("range_min", &SpatRaster::range_min)
        .def_property_readonly("range_max", &SpatRaster::range_max)

        .def("rowFromY",
            (std::vector<int64_t> (SpatRaster::*)(const std::vector<double>&))&SpatRaster::rowFromY)
        .def("colFromX",
            (std::vector<int64_t> (SpatRaster::*)(const std::vector<double>&))&SpatRaster::colFromX)
        .def("rowColFromCell",
            (std::vector<std::vector<int64_t>> (SpatRaster::*)(std::vector<double>))&SpatRaster::rowColFromCell)
        .def("cellFromXY",
            (std::vector<double> (SpatRaster::*)(std::vector<double>, std::vector<double>, double))&SpatRaster::cellFromXY)
        .def("cellFromRowCol",
            (std::vector<double> (SpatRaster::*)(std::vector<int64_t>, std::vector<int64_t>))&SpatRaster::cellFromRowCol)
        .def("yFromRow",
            (std::vector<double> (SpatRaster::*)(const std::vector<int64_t>&))&SpatRaster::yFromRow)
        .def("xFromCol",
            (std::vector<double> (SpatRaster::*)(const std::vector<int64_t>&))&SpatRaster::xFromCol)
        .def("xyFromCell",
            (std::vector<std::vector<double>> (SpatRaster::*)(std::vector<double>&))&SpatRaster::xyFromCell)
        .def("adjacent",      &SpatRaster::adjacent)

        .def("layerCor",      &SpatRaster::layerCor)
        .def("global_weighted_mean", &SpatRaster::global_weighted_mean)

        .def("subset",    &SpatRaster::subset)
        .def("selectRange",&SpatRaster::selRange)

        .def("setNames",  &SpatRaster::setNames)

        .def("arith_rast",
            (SpatRaster (SpatRaster::*)(SpatRaster, std::string, bool, SpatOptions&))&SpatRaster::arith)
        .def("arith_numb",
            (SpatRaster (SpatRaster::*)(std::vector<double>, std::string, bool, bool, SpatOptions&))&SpatRaster::arith)
        .def("arith_m",     &SpatRaster::arith_m)

        .def("rst_area",        &SpatRaster::rst_area)
        .def("sum_area",        &SpatRaster::sum_area)
        .def("sum_area_group",  &SpatRaster::sum_area_group)
        .def("surface_area",    &SpatRaster::surfaceArea)

        .def("as_points",       &SpatRaster::as_points)
        .def("cells_notna",     &SpatRaster::cells_notna)
        .def("cells_notna_novalues",&SpatRaster::cells_notna_novalues)
        .def("as_multipoints",  &SpatRaster::as_multipoints)
        .def("as_lines",        &SpatRaster::as_lines)
        .def("as_polygons",     &SpatRaster::as_polygons)
        .def("polygonize",      &SpatRaster::polygonize)

        .def("atan2",           &SpatRaster::atan_2)
        .def("bilinearValues",  &SpatRaster::bilinearValues)

        .def("patches",         &SpatRaster::clumps)
        .def("patches2",        &SpatRaster::patches)
        .def("collapse_sources",&SpatRaster::collapse_sources)
        .def("fill_range",      &SpatRaster::fill_range)
        .def("wmean_rast",
            (SpatRaster (SpatRaster::*)(SpatRaster, bool, SpatOptions&))(&SpatRaster::weighted_mean))
        .def("wmean_vect",
            (SpatRaster (SpatRaster::*)(std::vector<double>, bool, SpatOptions&))(&SpatRaster::weighted_mean))
        .def_property_readonly("origin",    &SpatRaster::origin)
        .def("boundaries",      &SpatRaster::edges)
        .def("buffer",          &SpatRaster::buffer)
        .def("gridDistance",    &SpatRaster::gridDistance)
        .def("costDistance",    &SpatRaster::costDistance)
        .def("rastDistance",    &SpatRaster::distance)
        .def("nearest",         &SpatRaster::nearest)
        .def("vectDistance",    &SpatRaster::distance_vector)
        .def("rastDirection",   &SpatRaster::direction)
        .def("vectDirectionRasterize",&SpatRaster::direction_rasterize)

        .def("get_tiles_ext",      &SpatRaster::get_tiles_extent)
        .def("get_tiles_ext_vect", &SpatRaster::get_tiles_extent_vect)
        .def("make_tiles",         &SpatRaster::make_tiles)
        .def("make_tiles_vect",    &SpatRaster::make_tiles_vect)
        .def("ext_from_rc",        &SpatRaster::ext_from_rc)

        .def("combineCats",  &SpatRaster::combineCats)
        .def("droplevels",   &SpatRaster::dropLevels)
        .def("clamp",        &SpatRaster::clamp)
        .def("clamp_raster", &SpatRaster::clamp_raster)
        .def("clamp_ts",     &SpatRaster::clamp_ts)
        .def("replaceValues",&SpatRaster::replaceValues)
        .def("lookup_classify", &SpatRaster::lookup_classify)
        .def("lookup_subst",    &SpatRaster::lookup_subst)
        .def("lookup_catalyze", &SpatRaster::lookup_catalyze)
        .def("classify",
            (SpatRaster (SpatRaster::*)(std::vector<double>, size_t, unsigned, bool, bool, double, bool, bool, bool, SpatOptions&))(&SpatRaster::reclassify))
        .def("selRange",     &SpatRaster::selRange)
        .def("separate",     &SpatRaster::separate)
        .def("sort",         &SpatRaster::sort)
        .def("intersect",    &SpatRaster::intersect)

        .def("cover",
            (SpatRaster (SpatRaster::*)(SpatRaster, std::vector<double>, SpatOptions&))(&SpatRaster::cover))
        .def("cover_self",
            (SpatRaster (SpatRaster::*)(std::vector<double>, SpatOptions&))(&SpatRaster::cover))

        .def("crop",         &SpatRaster::crop)
        .def("crop_mask",    &SpatRaster::cropmask)
        .def("cum",          &SpatRaster::cum)
        .def("disaggregate", &SpatRaster::disaggregate)
        .def("expand",       &SpatRaster::extend)
        .def("extractCell",  &SpatRaster::extractCell)
        .def("extractVector",&SpatRaster::extractVector)
        .def("extractVectorFlat",&SpatRaster::extractVectorFlat)
        .def("extractBuffer",&SpatRaster::extractBuffer)

        .def("flip",         &SpatRaster::flip)
        .def("focal",        &SpatRaster::focal)
        .def("focalValues",  &SpatRaster::focal_values)
        .def("count",        &SpatRaster::count)
        .def("freq",         &SpatRaster::freq)
        .def("geometry",     &SpatRaster::geometry)

        .def("get_aggregates",    &SpatRaster::get_aggregates)
        .def("get_aggregate_dims",&SpatRaster::get_aggregate_dims2)
        .def("globalTF",          &SpatRaster::globalTF)
        .def("mglobal",           &SpatRaster::mglobal)

        .def("initf",
            (SpatRaster (SpatRaster::*)(std::string, bool, SpatOptions&))(&SpatRaster::init))
        .def("initv",
            (SpatRaster (SpatRaster::*)(std::vector<double>, SpatOptions&))(&SpatRaster::init))
        .def("is_in",       &SpatRaster::is_in)
        .def("is_in_cells", &SpatRaster::is_in_cells)
        .def("countnan",    &SpatRaster::countnan)
        .def("is_wrapper",  &SpatRaster::is_wrapper)

        .def("logic_rast",
            (SpatRaster (SpatRaster::*)(SpatRaster, std::string, SpatOptions&))(&SpatRaster::logic))
        .def("logic_numb",
            (SpatRaster (SpatRaster::*)(std::vector<double>, std::string, SpatOptions&))(&SpatRaster::logic))
        .def("mask_self",
            (SpatRaster (SpatRaster::*)(SpatOptions&))(&SpatRaster::mask))
        .def("mask_raster",
            (SpatRaster (SpatRaster::*)(SpatRaster&, bool, std::vector<double>, double, SpatOptions&))(&SpatRaster::mask))
        .def("mask_vector",
            (SpatRaster (SpatRaster::*)(SpatVector&, bool, double, bool, SpatOptions&))(&SpatRaster::mask))
        .def("math",    &SpatRaster::math)
        .def("math2",   &SpatRaster::math2)
        .def("modal",   &SpatRaster::modal)
        .def("quantile",&SpatRaster::quantile)
        .def("rasterize",&SpatRaster::rasterize)
        .def("rasterizePointsV",
            (SpatRaster (SpatRaster::*)(SpatVector&, std::string, std::vector<double>&, bool, double, SpatOptions&))(&SpatRaster::rasterizePoints))
        .def("rasterizePointsXY",
            (SpatRaster (SpatRaster::*)(std::vector<double>&, std::vector<double>&, std::string, std::vector<double>&, bool, double, SpatOptions&))(&SpatRaster::rasterizePoints))
        .def("rasterizeLyr",    &SpatRaster::rasterizeLyr)
        .def("rasterizeGeom",   &SpatRaster::rasterizeGeom)
        .def("rasterizeWindow", &SpatRaster::rasterizeWindow)
        .def("wincircle",       &SpatRaster::win_circle)
        .def("winrect",         &SpatRaster::win_rect)
        .def("rgb2col",         &SpatRaster::rgb2col)
        .def("rgb2hsx",         &SpatRaster::rgb2hsx)
        .def("hsx2rgb",         &SpatRaster::hsx2rgb)

        .def("reverse",      &SpatRaster::reverse)
        .def("rotate",       &SpatRaster::rotate)

        .def("sampleStratifiedCells", &SpatRaster::sampleStratifiedCells)
        .def("sampleRegularRaster",   &SpatRaster::sampleRegularRaster)
        .def("sampleRowColRaster",    &SpatRaster::sampleRowColRaster)
        .def("sampleRegularValues",   &SpatRaster::sampleRegularValues)
        .def("sampleRowColValues",    &SpatRaster::sampleRowColValues)
        .def("sampleRowCol",          &SpatRaster::sampleRowCol)
        .def("sampleRandomRaster",    &SpatRaster::sampleRandomRaster)
        .def("sampleRandomValues",    &SpatRaster::sampleRandomValues)
        .def("scale",        &SpatRaster::scale)
        .def("scale_linear", &SpatRaster::scale_linear)
        .def("shift",        &SpatRaster::shift)
        .def("similarity",   &SpatRaster::similarity)
        .def("terrain",      &SpatRaster::terrain)
        .def("hillshade",    &SpatRaster::hillshade)
        .def("summary",      &SpatRaster::summary)
        .def("summary_numb", &SpatRaster::summary_numb)
        .def("transpose",    &SpatRaster::transpose)
        .def("trig",         &SpatRaster::trig)
        .def("trim1",        &SpatRaster::trim1)
        .def("trim",         &SpatRaster::trim2)
        .def("unique",       &SpatRaster::unique)
        .def("crosstab",     &SpatRaster::crosstab)
        .def("where",        &SpatRaster::where)
        .def("sieve",        &SpatRaster::sieveFilter)
        .def("view",         &SpatRaster::viewshed)
        .def("proximity",    &SpatRaster::proximity)
        .def("fillNA",       &SpatRaster::fillNA)
        .def("rectify",      &SpatRaster::rectify)
        .def("stretch",      &SpatRaster::stretch)
        .def("warp",         &SpatRaster::warper)
        .def("warp_by_util", &SpatRaster::warper_by_util)
        .def("resample",     &SpatRaster::resample)
        .def("zonal",        &SpatRaster::zonal)
        .def("zonal_weighted",&SpatRaster::zonal_weighted)
        .def("zonal_poly",   &SpatRaster::zonal_poly)
        .def("zonal_poly_table",   &SpatRaster::zonal_poly_table)
        .def("zonal_poly_weighted",&SpatRaster::zonal_poly_weighted)
        .def("watershed2",   &SpatRaster::watershed2)
        .def("pitfinder2",   &SpatRaster::pitfinder2)
        .def("NIDP2",        &SpatRaster::NIDP2)
        .def("flowAccu2",    &SpatRaster::flowAccu2)
        .def("flowAccu2_weight",&SpatRaster::flowAccu2_weight)
        .def("centroid",     &SpatRaster::centroid)
    ;


    // ── SpatRasterCollection ────────────────────────────────────────────────
    py::class_<SpatRasterCollection>(m, "SpatRasterCollection")
        .def(py::init<>())
        .def(py::init<std::string, std::vector<int>, bool,
                      std::vector<std::string>, bool, bool,
                      std::vector<std::string>>())
        .def_property("names",
            &SpatRasterCollection::get_names,
            &SpatRasterCollection::set_names)
        .def("deepcopy",    &SpatRasterCollection::deepCopy)
        .def("dims",        &SpatRasterCollection::dims)
        .def("extent",      &SpatRasterCollection::getExtent)
        .def("has_error",   &SpatRasterCollection::hasError)
        .def("has_warning", &SpatRasterCollection::hasWarning)
        .def("getError",    &SpatRasterCollection::getError)
        .def("getWarnings", &SpatRasterCollection::getWarnings)
        .def_readonly("x",  &SpatRasterCollection::ds)
        .def("length",      &SpatRasterCollection::size)
        .def("resize",      &SpatRasterCollection::resize)
        .def("erase",       &SpatRasterCollection::erase)
        .def("add",         &SpatRasterCollection::push_back)
        .def("merge",       &SpatRasterCollection::merge)
        .def("mosaic",      &SpatRasterCollection::mosaic)
        .def("morph",       &SpatRasterCollection::morph)
        .def("crop",        &SpatRasterCollection::crop)
        .def("addTag",      &SpatRasterCollection::addTag)
        .def("getTags",     &SpatRasterCollection::getTags)
        .def("make_vrt",    &SpatRasterCollection::make_vrt)
    ;


    // ── SpatRasterStack ─────────────────────────────────────────────────────
    py::class_<SpatRasterStack>(m, "SpatRasterStack")
        .def(py::init<>())
        .def(py::init<std::string, std::vector<int>, bool,
                      std::vector<std::string>, bool, bool,
                      std::vector<std::string>>())
        .def(py::init<SpatRaster, std::string, std::string, std::string>())
        .def("deepcopy",    &SpatRasterStack::deepCopy)
        .def("has_error",   &SpatRasterStack::hasError)
        .def("has_warning", &SpatRasterStack::hasWarning)
        .def("getError",    &SpatRasterStack::getError)
        .def("getWarnings", &SpatRasterStack::getWarnings)
        .def("readStart",   &SpatRasterStack::readStart)
        .def("readStop",    &SpatRasterStack::readStop)
        .def("readAll",     &SpatRasterStack::readAll)
        .def("nsds",        &SpatRasterStack::nsds)
        .def("ncol",        &SpatRasterStack::ncol)
        .def("nrow",        &SpatRasterStack::nrow)
        .def("nlyr",        &SpatRasterStack::nlyr)
        .def("res",         &SpatRasterStack::resolution)
        .def("ext",         &SpatRasterStack::getExtent)
        .def("filenames",   &SpatRasterStack::filenames)
        .def("get_crs",     &SpatRasterStack::getSRS)
        .def_property("names",
            &SpatRasterStack::get_names, &SpatRasterStack::set_names)
        .def_property("long_names",
            &SpatRasterStack::get_longnames, &SpatRasterStack::set_longnames)
        .def_property("units",
            &SpatRasterStack::get_units, &SpatRasterStack::set_units)
        .def("add",       &SpatRasterStack::push_back)
        .def("resize",    &SpatRasterStack::resize)
        .def("summary",   &SpatRasterStack::summary)
        .def("summary_numb",&SpatRasterStack::summary_numb)
        .def("getsds",    &SpatRasterStack::getsds)
        .def("replace",   &SpatRasterStack::replace)
        .def("subset",    &SpatRasterStack::subset)
        .def("collapse",  &SpatRasterStack::collapse)
        .def("extractCell",  &SpatRasterStack::extractCell)
        .def("extractVector",&SpatRasterStack::extractVector)
        .def("crop",      &SpatRasterStack::crop)
        .def("addTag",    &SpatRasterStack::addTag)
        .def("getTags",   &SpatRasterStack::getTags)
    ;
}
