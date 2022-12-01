#include <Rcpp.h>
//#include "spatRaster.h"
#include "spatRasterMultiple.h"
#include "spatGraph.h"
#include <memory> //std::addressof
#include "NA.h"
#include "spatTime.h"
//#include "spatVector2.h"

//static void SpatRaster_finalizer( SpatRaster* ptr ){
//}

Rcpp::List getBlockSizeR(SpatRaster* r, unsigned n, double frac) {
	SpatOptions opt;
	opt.ncopies = n;
	opt.set_memfrac(frac);
    BlockSize bs = r->getBlockSize(opt);
	Rcpp::List L = Rcpp::List::create(Rcpp::Named("row") = bs.row, Rcpp::Named("nrows") = bs.nrows, Rcpp::Named("n") = bs.n);
	return(L);
}


Rcpp::List getBlockSizeWrite(SpatRaster* r) {
    BlockSize bs = r->bs;
	Rcpp::List L = Rcpp::List::create(Rcpp::Named("row") = bs.row, Rcpp::Named("nrows") = bs.nrows, Rcpp::Named("n") = bs.n);
	return(L);
}


bool sameObject(SpatRaster* a, SpatRaster* b) {
	return a == b;
}



Rcpp::List getDataFrame(SpatDataFrame* v) {
	unsigned n = v->ncol();
	Rcpp::List out(n);
	if (n == 0) {
		return(out);
	}

	std::vector<std::string> nms = v->names;
	std::vector<unsigned> itype = v->itype;
	for (size_t i=0; i < n; i++) {
		if (itype[i] == 0) {
			out[i] = v->getD(i);
		} else if (itype[i] == 1) {
			Rcpp::NumericVector iv = Rcpp::wrap(v->getI(i));
			long longNA = NA<long>::value;
			for (R_xlen_t j=0; j<iv.size(); j++) {
				if (iv[j] == longNA) {
					iv[j] = NA_REAL;
				}
			}
			out[i] = iv;

		} else if (itype[i] == 2){
			Rcpp::CharacterVector s = Rcpp::wrap(v->getS(i));
			for (R_xlen_t j=0; j<s.size(); j++) {
				if (s[j] == "____NA_+") {
					s[j] = NA_STRING;
				}
			}
			out[i] = s;
		} else if (itype[i] == 3){
			std::vector<int8_t> b = v->getB(i);
			Rcpp::NumericVector d(b.size());
			for (size_t j=0; j<b.size(); j++) {
				if (b[j] > 1) {
					d[j] = NA_REAL;
				} else {
					d[j] = b[j];
				}
			}
			out[i] = d;
		} else if (itype[i] == 4){
			SpatTime_v tx = v->getT(i);
			Rcpp::NumericVector tv = Rcpp::wrap(tx.x);
			SpatTime_t timeNA = NA<SpatTime_t>::value;
			for (R_xlen_t j=0; j<tv.size(); j++) {
				if (tv[j] == timeNA) {
					tv[j] = NA_REAL;
				}
			}
			out[i] = tv;
		} else if (itype[i] == 5){
			out[i] = v->getF(i);
		}
	}
	out.names() = nms;
	// todo: deal with NAs in int and str
	return out;
//  Rcpp::df is nice, but no of variables is <= 20,
//  and no "stringsAsFactors"=false
//	Rcpp::DataFrame result(out);
//	result.attr("names") = v->names();
//	return result;
}


Rcpp::List getVectorAttributes(SpatVector* v) {
	SpatDataFrame df = v->df;
	Rcpp::List lst = getDataFrame(&df);
	return lst;
}


/*
Rcpp::List getRasterAttributes(SpatRaster* x) {
	Rcpp::List lst;
	if (x->nlyr() > 0) {
		SpatDataFrame df = x->source[0].cats[0];
		lst = getDataFrame(&df);
	}
	return lst;
}
*/


Rcpp::DataFrame get_geometryDF(SpatVector* v) {
	SpatDataFrame df = v->getGeometryDF();

	Rcpp::DataFrame out = Rcpp::DataFrame::create(
			Rcpp::Named("id") = df.iv[0],
			Rcpp::Named("part") = df.iv[1],
			Rcpp::Named("x") = df.dv[0],
			Rcpp::Named("y") = df.dv[1],
			Rcpp::Named("hole") = df.iv[2]
	);
	return out;
}

std::vector<std::vector<std::vector<Rcpp::DataFrame>>> get_geometryList(SpatVector* v, const std::string xnm, const std::string ynm) {
	size_t ni = v->nrow();
	std::vector<std::vector<std::vector<Rcpp::DataFrame>>> out(ni);
	for (size_t i=0; i < ni; i++) {
		SpatGeom g = v->getGeom(i);
		size_t nj = g.size();
		if (nj == 0) { // empty
			continue;
		}
		out[i].resize(nj);
		for (size_t j=0; j < nj; j++) {
			SpatPart p = g.getPart(j);
			size_t nk = p.nHoles();
			out[i][j].reserve(nk + 1);
			Rcpp::DataFrame m = Rcpp::DataFrame::create(
				Rcpp::Named(xnm) = p.x,
				Rcpp::Named(ynm) = p.y
			);
			out[i][j].push_back(m);
			for (size_t k=0; k<nk ; k++) {
				SpatHole h = p.getHole(k);
				Rcpp::DataFrame m = Rcpp::DataFrame::create(
					Rcpp::Named(xnm) = h.x,
					Rcpp::Named(ynm) = h.y
				);
				out[i][j].push_back(m);
			}
		}
	}
	return out;
}


RCPP_EXPOSED_CLASS(SpatFactor)
RCPP_EXPOSED_CLASS(SpatTime_v)
RCPP_EXPOSED_CLASS(SpatSRS)
RCPP_EXPOSED_CLASS(SpatExtent)
RCPP_EXPOSED_CLASS(SpatMessages)
RCPP_EXPOSED_CLASS(SpatOptions)
RCPP_EXPOSED_CLASS(SpatCategories)
RCPP_EXPOSED_CLASS(SpatDataFrame)
RCPP_EXPOSED_CLASS(SpatWindow)
RCPP_EXPOSED_CLASS(SpatRasterSource)
RCPP_EXPOSED_CLASS(SpatRaster)
RCPP_EXPOSED_CLASS(SpatRasterCollection)
RCPP_EXPOSED_CLASS(SpatRasterStack)
RCPP_EXPOSED_CLASS(SpatVector)
RCPP_EXPOSED_CLASS(SpatVectorProxy)
RCPP_EXPOSED_CLASS(SpatVectorCollection)
RCPP_EXPOSED_CLASS(SpatGraph)
//RCPP_EXPOSED_CLASS(SpatVector2)

RCPP_MODULE(spat){

    using namespace Rcpp;

    class_<SpatTime_v>("SpatTime_v")
		.constructor()
		.field("step", &SpatTime_v::step)
		.field("zone", &SpatTime_v::zone)
		.field("x", &SpatTime_v::x)
	;

    class_<SpatFactor>("SpatFactor")
		.constructor()
		.constructor<std::vector<unsigned>, std::vector<std::string>>()
		.field("values", &SpatFactor::v)
		.field("labels", &SpatFactor::labels)
	;


    class_<SpatSRS>("SpatSRS")
		.constructor()
		.method("set", &SpatSRS::set)
		.method("is_lonlat", &SpatSRS::is_lonlat)
		.method("to_meter", &SpatSRS::to_meter)
	;

    class_<SpatGraph>("SpatGraph")
		.constructor()
	;

    class_<SpatExtent>("SpatExtent")
		.constructor()
		.constructor<double, double, double, double>()
		.method("deepcopy", &SpatExtent::deepCopy, "deepCopy")
		.property("vector", &SpatExtent::asVector)
		.property("valid", &SpatExtent::valid)
		.method("align", &SpatExtent::align, "align")
		.method("intersect", &SpatExtent::intersect, "intersect")
		.method("as.points", &SpatExtent::asPoints, "as.points")
		.method("ceil",  &SpatExtent::ceil,  "ceil")
		.method("compare", &SpatExtent::compare, "compare")
		.method("floor", &SpatExtent::floor, "floor")
		.method("round", &SpatExtent::round, "round")
		.method("union", &SpatExtent::unite, "union")
		.method("sampleRandom", &SpatExtent::sampleRandom)
		.method("sampleRegular", &SpatExtent::sampleRegular)
		.method("sample", &SpatExtent::test_sample)
	;

/*
    class_<SpatWindow>("SpatWindow")
		.field_readonly("full_extent", &SpatWindow::full_extent)
		.field_readonly("full_nrow", &SpatWindow::full_nrow)
		.field_readonly("full_ncol", &SpatWindow::full_ncol)
		.field_readonly("off_row", &SpatWindow::off_row)
		.field_readonly("off_col", &SpatWindow::off_col)
		.field_readonly("expand", &SpatWindow::expand)
	;
*/

    class_<SpatMessages>("SpatMessages")
		.constructor()
		//.field("success", &SpatMessages::success)
		.field("has_error", &SpatMessages::has_error)
		.field("has_warning", &SpatMessages::has_warning)
		.method("getError", &SpatMessages::getError)
		.method("getWarnings", &SpatMessages::getWarnings)
	;

    class_<SpatOptions>("SpatOptions")
		.constructor()
		.method("deepcopy", &SpatOptions::deepCopy, "deepCopy")
		.property("tempdir", &SpatOptions::get_tempdir, &SpatOptions::set_tempdir )
		.property("memfrac", &SpatOptions::get_memfrac, &SpatOptions::set_memfrac )
		.property("memmax", &SpatOptions::get_memmax, &SpatOptions::set_memmax )
		.property("memmin", &SpatOptions::get_memmin, &SpatOptions::set_memmin )
		.property("tolerance", &SpatOptions::get_tolerance, &SpatOptions::set_tolerance )
		.property("filenames", &SpatOptions::get_filenames, &SpatOptions::set_filenames )
		.property("filetype", &SpatOptions::get_filetype, &SpatOptions::set_filetype )
		.property("datatype", &SpatOptions::get_datatype, &SpatOptions::set_datatype )
		.property("verbose", &SpatOptions::get_verbose, &SpatOptions::set_verbose )
		.property("NAflag", &SpatOptions::get_NAflag, &SpatOptions::set_NAflag )
		//.property("ncdfcopy", &SpatOptions::get_ncdfcopy, &SpatOptions::set_ncdfcopy )
		.property("statistics", &SpatOptions::get_statistics, &SpatOptions::set_statistics )
		.property("overwrite", &SpatOptions::get_overwrite, &SpatOptions::set_overwrite )
		//.property("append", &SpatOptions::get_append, &SpatOptions::set_append )
		.field("datatype_set", &SpatOptions::datatype_set)
		.field("threads", &SpatOptions::threads)
		.property("progress", &SpatOptions::get_progress, &SpatOptions::set_progress)
		.field("progressbar", &SpatOptions::progressbar)
		.property("ncopies", &SpatOptions::get_ncopies, &SpatOptions::set_ncopies)

		.property("def_filetype", &SpatOptions::get_def_filetype, &SpatOptions::set_def_filetype )
		.property("def_datatype", &SpatOptions::get_def_datatype, &SpatOptions::set_def_datatype )
		.property("def_verbose", &SpatOptions::get_def_verbose, &SpatOptions::set_def_verbose )

		//.property("def_bandorder", &SpatOptions::get_def_bandorder, &SpatOptions::set_def_bandorder )

		.property("todisk", &SpatOptions::get_todisk, &SpatOptions::set_todisk)
		.field("messages", &SpatOptions::msg, "messages")
		.field("gdal_options", &SpatOptions::gdal_options)
		.field("pid", &SpatOptions::pid )
		.field("names", &SpatOptions::names)
		.property("steps", &SpatOptions::get_steps, &SpatOptions::set_steps)
	//	.property("overwrite", &SpatOptions::set_overwrite, &SpatOptions::get_overwrite )
		//.field("gdaloptions", &SpatOptions::gdaloptions)

		.property("scale", &SpatOptions::get_scale, &SpatOptions::set_scale)
		.property("offset", &SpatOptions::get_offset, &SpatOptions::set_offset)

	;

    class_<SpatDataFrame>("SpatDataFrame")
		.constructor()

		.field_readonly("itype", &SpatDataFrame::itype)
		.field_readonly("iplace", &SpatDataFrame::iplace)

		.property("names", &SpatDataFrame::get_names, &SpatDataFrame::set_names)
		.property("nrow", &SpatDataFrame::nrow, &SpatDataFrame::resize_rows, "nrow")
		.property("ncol", &SpatDataFrame::ncol, &SpatDataFrame::resize_cols, "ncol")

		.method("has_error", &SpatDataFrame::hasError)
		.method("has_warning", &SpatDataFrame::hasWarning)
		.method("getWarnings", &SpatDataFrame::getWarnings)
		.method("getError", &SpatDataFrame::getError)

		.method("add_column_double", (bool (SpatDataFrame::*)(std::vector<double>, std::string name))( &SpatDataFrame::add_column))
		.method("add_column_long", (bool (SpatDataFrame::*)(std::vector<long>, std::string name))( &SpatDataFrame::add_column))
		.method("add_column_string", (bool (SpatDataFrame::*)(std::vector<std::string>, std::string name))( &SpatDataFrame::add_column))
		.method("add_column_factor", (bool (SpatDataFrame::*)(SpatFactor, std::string name))( &SpatDataFrame::add_column))
		.method("add_column_bool", &SpatDataFrame::add_column_bool)
		.method("add_column_time", &SpatDataFrame::add_column_time)

		.method("remove_column", (bool (SpatDataFrame::*)(std::string field))( &SpatDataFrame::remove_column))
		.method("remove_column", (bool (SpatDataFrame::*)(int i))( &SpatDataFrame::remove_column))
		.method("get_datatypes", &SpatDataFrame::get_datatypes)
		.method("get_timezones", &SpatDataFrame::get_timezones)
		.method("get_timesteps", &SpatDataFrame::get_timesteps)

		.method("subset_rows", (SpatDataFrame (SpatDataFrame::*)(std::vector<unsigned>))( &SpatDataFrame::subset_rows))
		.method("subset_cols", (SpatDataFrame (SpatDataFrame::*)(std::vector<unsigned>))( &SpatDataFrame::subset_cols))

		.method("remove_rows", &SpatDataFrame::remove_rows)

		.method("cbind", &SpatDataFrame::cbind)
		.method("rbind", &SpatDataFrame::rbind)
		.method("values", &getDataFrame)
		.method("unique", &SpatDataFrame::unique)
		//.method("one_string", &SpatDataFrame::one_string)

		.method("write", &SpatDataFrame::write_dbf)
		.field("messages", &SpatDataFrame::msg)
		.method("strwidth", &SpatDataFrame::strwidth)
	;


    class_<SpatVectorCollection>("SpatVectorCollection")
		.constructor()

		//.property("names", &SpatVectorCollection::get_names, &SpatVectorCollection::set_names)
		.method("deepcopy", &SpatVectorCollection::deepCopy, "deepCopy")

		.method("size", &SpatVectorCollection::size)
		.method("get", &SpatVectorCollection::get)
		.method("push_back", &SpatVectorCollection::push_back)
		.method("subset", &SpatVectorCollection::subset)
		.method("replace", &SpatVectorCollection::replace)
		.method("append", &SpatVectorCollection::append)

		.method("has_error", &SpatVectorCollection::hasError)
		.method("has_warning", &SpatVectorCollection::hasWarning)
		.method("getWarnings", &SpatVectorCollection::getWarnings)
		.method("getError", &SpatVectorCollection::getError)
		.method("from_hex_col", &SpatVectorCollection::from_hex_col)
		.method("setNames", &SpatVectorCollection::setNames, "setNames" )
		.property("names", &SpatVectorCollection::getNames)

	;


    class_<SpatCategories>("SpatCategories")
		.constructor()
		.field_readonly("df", &SpatCategories::d)
		.field("index", &SpatCategories::index)
		.method("combine", &SpatCategories::combine)
	;


    class_<SpatVector>("SpatVector")
		.constructor()
		.constructor<SpatExtent, std::string>()
		.constructor<std::vector<std::string>>()
	;


    class_<SpatVector>("SpatVector")
		.constructor()
		.constructor<SpatExtent, std::string>()
		.constructor<std::vector<std::string>>()

		.method("deepcopy", &SpatVector::deepCopy)
		.method("wkt", &SpatVector::wkt)
		.method("wkb", &SpatVector::wkb)
		.method("hex", &SpatVector::hex)
		.method("from_hex", &SpatVector::from_hex)
		.method("make_nodes", &SpatVector::make_nodes)
		.method("boundary", &SpatVector::boundary)
		.method("polygonize", &SpatVector::polygonize)
		.method("normalize", &SpatVector::normalize)
		.method("normalize_longitude", &SpatVector::normalize_longitude)
		.method("rotate_longitude", &SpatVector::rotate_longitude)
		.method("line_merge", &SpatVector::line_merge)
		.method("simplify", &SpatVector::simplify)
		.method("thin", &SpatVector::thin)
		//.method("shared_paths", &SpatVector::shared_paths)
		.method("shared_paths", (SpatVector (SpatVector::*)())( &SpatVector::shared_paths))
		.method("shared_paths2", (SpatVector (SpatVector::*)(SpatVector))( &SpatVector::shared_paths))
		.method("snap", &SpatVector::snap)
		.method("snapto", &SpatVector::snapto)
		.method("spatial_index_2d", &SpatVector::index_2d)
		.method("spatial_index_sparse", &SpatVector::index_sparse)

		.field_readonly("is_proxy", &SpatVector::is_proxy )
		.field_readonly("read_query", &SpatVector::read_query )
		.field_readonly("read_extent", &SpatVector::read_extent )
		.field_readonly("geom_count", &SpatVector::geom_count)
		.field_readonly("source", &SpatVector::source)
		.field_readonly("layer", &SpatVector::source_layer)
		.field_readonly("df", &SpatVector::df )

		.method("has_error", &SpatVector::hasError)
		.method("has_warning", &SpatVector::hasWarning)
		.method("getError", &SpatVector::getError)
		.method("getWarnings", &SpatVector::getWarnings)

		.method("coordinates", &SpatVector::coordinates)
		.method("get_geometry", &SpatVector::getGeometry)
		.method("get_geometryDF", &get_geometryDF)
		.method("get_geometryList", &get_geometryList)
		.method("linesList", &SpatVector::linesList)
		.method("polygonsList", &SpatVector::polygonsList)
		.method("linesNA", &SpatVector::linesNA)
		.method("linedistlonlat", &SpatVector::linedistLonLat)

		.method("add_column_empty", (void (SpatVector::*)(unsigned dtype, std::string name))( &SpatVector::add_column))
		.method("add_column_double", (bool (SpatVector::*)(std::vector<double>, std::string name))( &SpatVector::add_column))
		.method("add_column_long", (bool (SpatVector::*)(std::vector<long>, std::string name))( &SpatVector::add_column))
		.method("add_column_string", (bool (SpatVector::*)(std::vector<std::string>, std::string name))( &SpatVector::add_column))
		.method("add_column_factor", &SpatVector::add_column_factor)
		.method("add_column_bool", &SpatVector::add_column_bool)
		.method("add_column_time", &SpatVector::add_column_time)
		.method("remove_column", (bool (SpatVector::*)(std::string field))( &SpatVector::remove_column))
		.method("remove_column", (bool (SpatVector::*)(int i))( &SpatVector::remove_column))
		.method("remove_df", &SpatVector::remove_df)
		.method("set_df", &SpatVector::set_df)
		.method("get_datatypes", &SpatVector::get_datatypes)

		.method("set_holes", &SpatVector::set_holes)
		.method("get_holes", &SpatVector::get_holes)
		.method("remove_holes", &SpatVector::remove_holes)
		.method("append", &SpatVector::append)
		.method("cbind", &SpatVector::cbind)
		.method("rbind", &SpatVector::append)

		.method("area", &SpatVector::area)
		.method("as_lines", &SpatVector::as_lines)
		.method("as_points", &SpatVector::as_points)
		.method("couldBeLonLat", &SpatVector::could_be_lonlat)
		.method("get_crs", &SpatVector::getSRS)

		.method("set_crs", (bool (SpatVector::*)(std::string crs))( &SpatVector::setSRS))
		//.method("prj", &SpatVector::getPRJ)

		.method("distance_self", (std::vector<double> (SpatVector::*)(bool, std::string))( &SpatVector::distance))
		.method("distance_other", (std::vector<double> (SpatVector::*)(SpatVector, bool, std::string))( &SpatVector::distance))

		.method("point_distance", &SpatVector::pointdistance)

		.method("geosdist_self", (std::vector<double> (SpatVector::*)(bool, std::string))( &SpatVector::geos_distance))
		.method("geosdist_other", (std::vector<double> (SpatVector::*)(SpatVector, bool, std::string))( &SpatVector::geos_distance))

		.method("extent", &SpatVector::getExtent)
		.method("getDF", &getVectorAttributes)
		.method("getGeometryWKT", &SpatVector::getGeometryWKT)
		.method("isLonLat", &SpatVector::is_lonlat)
		.method("length", &SpatVector::length)
//		.field("srs", &SpatVector::srs, "srs")
		.field("messages", &SpatVector::msg)
		.property("names", &SpatVector::get_names, &SpatVector::set_names)
		.method("nrow", &SpatVector::nrow)
		.method("ncol", &SpatVector::ncol)
		.method("project", &SpatVector::project)
		.method("project_xy", &SpatVector::project_xy)
		.method("read", &SpatVector::read)
		.method("setGeometry", &SpatVector::setGeometry)
		.method("setPointsXY", &SpatVector::setPointsGeometry)
		.method("setPointsDF", &SpatVector::setPointsDF)
		.method("size", &SpatVector::size)
		.method("subset_cols", ( SpatVector (SpatVector::*)(std::vector<int>))( &SpatVector::subset_cols ))
		.method("subset_rows", ( SpatVector (SpatVector::*)(std::vector<int>))( &SpatVector::subset_rows ))
		.method("remove_rows", &SpatVector::remove_rows)
		.method("type", &SpatVector::type)
		.method("nullGeoms", &SpatVector::nullGeoms)

		.method("write", &SpatVector::write)
		.method("delete_layers", &SpatVector::delete_layers)
		.method("layer_names", &SpatVector::layer_names)

//		.method("bienvenue", &SpatVector::bienvenue)
//		.method("allerretour", &SpatVector::allerretour)
		.method("geos_isvalid", &SpatVector::geos_isvalid)
		.method("geos_isvalid_msg", &SpatVector::geos_isvalid_msg)

		.method("aggregate", ( SpatVector (SpatVector::*)(std::string, bool))( &SpatVector::aggregate ))
		.method("aggregate_nofield", ( SpatVector (SpatVector::*)(bool))( &SpatVector::aggregate ))

		.method("disaggregate", &SpatVector::disaggregate)
		.method("buffer", &SpatVector::buffer)
		.method("centroid", &SpatVector::centroid)
		.method("point_on_surface", &SpatVector::point_on_surface)
		.method("make_valid2", &SpatVector::make_valid2)
		.method("flip", &SpatVector::flip)
		.method("transpose", &SpatVector::transpose)
		.method("shift", &SpatVector::shift)
		.method("rescale", &SpatVector::rescale)
		.method("rotate", &SpatVector::rotate)
		.method("erase_agg", &SpatVector::erase_agg)
		.method("erase", ( SpatVector (SpatVector::*)(SpatVector))( &SpatVector::erase ))
		.method("erase_self", ( SpatVector (SpatVector::*)(bool))( &SpatVector::erase ))
		.method("elongate", &SpatVector::elongate)
		.method("gaps", &SpatVector::gaps)

		.method("symdif", &SpatVector::symdif)
		.method("cover", &SpatVector::cover)
		.method("union", ( SpatVector (SpatVector::*)(SpatVector))( &SpatVector::unite ))
		.method("union_self", ( SpatVector (SpatVector::*)())( &SpatVector::unite ))
		.method("union_unary", &SpatVector::unaryunion)
		.method("intersect", &SpatVector::intersect)
		.method("delaunay", &SpatVector::delaunay)
		.method("voronoi", &SpatVector::voronoi)
		.method("hull", &SpatVector::hull)

		.method("width", &SpatVector::width)
		.method("clearance", &SpatVector::clearance)
		.method("mask", &SpatVector::mask)

		.method("is_related", &SpatVector::is_related)

		.method("related_between", ( std::vector<std::vector<double>> (SpatVector::*)(SpatVector, std::string, bool))( &SpatVector::which_relate))
		.method("related_within", ( std::vector<std::vector<double>> (SpatVector::*)(std::string, bool))( &SpatVector::which_relate))

//		.method("relate_first", &SpatVector::relateFirst)
		.method("relate_between", ( std::vector<int> (SpatVector::*)(SpatVector, std::string, bool, bool))( &SpatVector::relate ))
		.method("relate_within", ( std::vector<int> (SpatVector::*)(std::string, bool))( &SpatVector::relate ))
		.method("equals_between", ( std::vector<unsigned> (SpatVector::*)(SpatVector, double))( &SpatVector::equals_exact ))
		.method("equals_within", ( std::vector<unsigned> (SpatVector::*)(bool, double))( &SpatVector::equals_exact ))


		.method("crop_ext", ( SpatVector (SpatVector::*)(SpatExtent))( &SpatVector::crop ))
		.method("crop_vct", ( SpatVector (SpatVector::*)(SpatVector))( &SpatVector::crop ))

		.method("near_between", (SpatVector (SpatVector::*)(SpatVector, bool))( &SpatVector::nearest_point))
		.method("near_within", (SpatVector (SpatVector::*)())( &SpatVector::nearest_point))
		.method("near_geom", &SpatVector::nearest_geometry)

		//.method("knearest", &SpatVector::knearest)

		.method("split", &SpatVector::split)

		.method("sample", &SpatVector::sample)
		.method("sampleGeom", &SpatVector::sample_geom)
		.method("remove_duplicate_nodes", &SpatVector::remove_duplicate_nodes)

		.method("cross_dateline", &SpatVector::cross_dateline)
		.method("fix_lonlat_overflow", &SpatVector::fix_lonlat_overflow)

		.method("densify", &SpatVector::densify)
		.method("round", &SpatVector::round)
	;


//    class_<SpatRasterSource>("SpatRasterSource")
//		.field_readonly("cats", &SpatRasterSource::cats)
///		.field_readonly("has_scale_offset", &SpatRasterSource::has_scale_offset)
//		.field_readonly("scale", &SpatRasterSource::scale)
//		.field_readonly("offset", &SpatRasterSource::offset)

//		.field_readonly("time", &SpatRasterSource::time)
//		.field("srs", &SpatRasterSource::srs, "srs")

		//.field_readonly("memory", &SpatRasterSource::memory)
	//	.field_readonly("filename", &SpatRasterSource::filename)
		//.field_readonly("driver", &SpatRasterSource::driver)
		//.field_readonly("nrow", &SpatRasterSource::nrow)
		//.field_readonly("ncol", &SpatRasterSource::ncol)
		//.field_readonly("nlyr", &SpatRasterSource::nlyr)
//		.field_readonly("extent", &SpatRasterSource::extent)
//		.field_readonly("hasWindow", &SpatRasterSource::hasWindow)
//		.field_readonly("window", &SpatRasterSource::window)
		//.field_readonly("layers", &SpatRasterSource::layers)
		//.field_readonly("nlyrfile", &SpatRasterSource::nlyrfile)
		//.field_readonly("flipped", &SpatRasterSource::flipped)
		//.field_readonly("rotated", &SpatRasterSource::rotated)
//		.field_readonly("parameters_changed", &SpatRasterSource::parameters_changed)
//;


    class_<SpatVectorProxy>("SpatVectorProxy")
		.constructor()
		.field("v", &SpatVectorProxy::v )
		.method("deepcopy", &SpatVectorProxy::deepCopy, "deepCopy")
	;



    class_<SpatRaster>("SpatRaster")
		.constructor()
	 // .constructor<std::string, int>()
	    .constructor<std::vector<std::string>, std::vector<int>, std::vector<std::string>, bool, std::vector<std::string>, std::vector<std::string>, std::vector<size_t>>()
		.constructor<std::vector<unsigned>, std::vector<double>, std::string>()
		//.finalizer(&SpatRaster_finalizer)

		//.method("fromFiles", &SpatRaster::fromFiles)
		.method("has_error", &SpatRaster::hasError)
		.method("has_warning", &SpatRaster::hasWarning)
		.method("getError", &SpatRaster::getError)
		.method("getWarnings", &SpatRaster::getWarnings)
		.method("getMessage", &SpatRaster::getMessage)

		//.field("name", &SpatRaster::name)
		.method("getFileBlocksize", &SpatRaster::getFileBlocksize)

		.method("sources_to_disk", &SpatRaster::sources_to_disk, "sources_to_disk")
		.method("mem_needs", &SpatRaster::mem_needs, "mem_needs")
		.method("spatinit", &SpatRaster::gdalogrproj_init, "init")
		.method("addSource", &SpatRaster::addSource, "addSource")
		.method("replace", &SpatRaster::replace, "replace")
		.method("combineSources", &SpatRaster::combineSources, "combineSources")
		.method("compare_geom", &SpatRaster::compare_geom, "compare_geom")
		.method("couldBeLonLat", &SpatRaster::could_be_lonlat, "couldBeLonLat")
		.method("deepcopy", &SpatRaster::deepCopy, "deepCopy")
		.method("hardcopy", &SpatRaster::hardCopy)
		.method("get_crs", &SpatRaster::getSRS)
		.method("set_crs", (bool (SpatRaster::*)(std::string crs))( &SpatRaster::setSRS))
		//.field_readonly("prj", &SpatRaster::prj)
		.property("extent", &SpatRaster::getExtent, &SpatRaster::setExtent )

		.method("setWindow", &SpatRaster::setWindow, "")
		.method("removeWindow", &SpatRaster::removeWindow, "")
		.method("hasWindow", &SpatRaster::hasWindow, "")

		//.method("getRasterAtt", &getRasterAttributes, "get attributes")

		.method("filenames", &SpatRaster::filenames )

		.field_readonly("rgb", &SpatRaster::rgb)
		.field_readonly("rgbtype", &SpatRaster::rgbtype)
		.method("setRGB", &SpatRaster::setRGB)
		.method("removeRGB", &SpatRaster::removeRGB)
		.method("getRGB", &SpatRaster::getRGB)

		//.method("setAttIndex", &SpatRaster::setAttrIndex, "setAttrIndex")
		//.method("getAttIndex", &SpatRaster::getAttrIndex, "getAttrIndex")
		//.method("hasAttributes", &SpatRaster::hasAttributes, "hasAttributes")
		//.method("getAttributes", &SpatRaster::getAttributes, "getAttributes")
		//.method("setAttributes", &SpatRaster::setAttributes, "setAttributes")
		//.method("createAttributes", &SpatRaster::createAttributes, "createAttributes")

		.method("makeCategorical", &SpatRaster::makeCategorical)
		.method("createCategories", &SpatRaster::createCategories)
		.method("hasCategories", &SpatRaster::hasCategories)
		.method("getCategories", &SpatRaster::getCategories)
		.method("setCategories", &SpatRaster::setCategories)
		.method("removeCategories", &SpatRaster::removeCategories)
		.method("getLabels", &SpatRaster::getLabels)
		.method("setLabels", &SpatRaster::setLabels)
		.method("getCatIndex", &SpatRaster::getCatIndex)
		.method("setCatIndex", &SpatRaster::setCatIndex)
		.method("getScaleOffset", &SpatRaster::getScaleOffset)
		.method("setScaleOffset", &SpatRaster::setScaleOffset)

		.method("hasColors", &SpatRaster::hasColors)
		.method("getColors", &SpatRaster::getColors)
		.method("setColors", &SpatRaster::setColors)
		.method("removeColors", &SpatRaster::removeColors)

		// need to keep this for raster
		.property("dataType", &SpatRaster::dataType)
		.method("getDataType", &SpatRaster::getDataType)
		.method("valueType", &SpatRaster::getValueType)
		.method("setValueType", &SpatRaster::setValueType)
		.property("hasRange", &SpatRaster::hasRange )
		.property("hasValues", &SpatRaster::hasValues )
		.property("inMemory", &SpatRaster::inMemory )
		.method("isLonLat", &SpatRaster::is_lonlat, "isLonLat")
		.method("isGlobalLonLat", &SpatRaster::is_global_lonlat, "isGlobalLonLat")

		.property("names", &SpatRaster::getNames)
		.method("get_sourcenames_long", &SpatRaster::getLongSourceNames)
		.method("set_sourcenames_long", &SpatRaster::setLongSourceNames)
		.method("get_sourcenames", &SpatRaster::getSourceNames)
		.method("set_sourcenames", &SpatRaster::setSourceNames)

		.method("setNAflag", &SpatRaster::setNAflag)
		.method("getNAflag", &SpatRaster::getNAflag)

		.property("hasUnit", &SpatRaster::hasUnit)
		.property("hasTime", &SpatRaster::hasTime)
		.property("time", &SpatRaster::getTime)
		.property("timestep", &SpatRaster::getTimeStep)
		.property("timezone", &SpatRaster::getTimeZone)
		.method("settime", &SpatRaster::setTime)
		//.property("timestr", &SpatRaster::getTimeStr)

		.property("depth", &SpatRaster::getDepth)
		.method("set_depth", &SpatRaster::setDepth)
		.property("units", &SpatRaster::getUnit)
		.method("set_units", &SpatRaster::setUnit)

		.method("size", &SpatRaster::size, "size")
		.method("nrow", &SpatRaster::nrow, "nrow")
		.method("ncol", &SpatRaster::ncol, "ncol")
		.method("nsrc", &SpatRaster::nsrc, "nsrc" )
		.field("messages", &SpatRaster::msg, "messages")
		.method("nlyrBySource", &SpatRaster::nlyrBySource, "nlyrBySource" )
		.method("lyrsBySource", &SpatRaster::lyrsBySource, "lyrsBySource" )
		.method("getBands", &SpatRaster::getBands, "getBands" )

		.method("nlyr", &SpatRaster::nlyr, "nlyr" )
		.property("origin", &SpatRaster::origin)
		.property("range_min", &SpatRaster::range_min )
		.property("range_max", &SpatRaster::range_max )
		.property("res", &SpatRaster::resolution)

// only if SpatRasterSource is exposed
//		.field_readonly("source", &SpatRaster::source )

		.method("collapse_sources", &SpatRaster::collapse_sources)
		.method("make_vrt", &SpatRaster::make_vrt)

		.method("dense_extent", &SpatRaster::dense_extent)
		.method("setNames", &SpatRaster::setNames)
		.method("setTime", &SpatRaster::setTime)
		.method("setDepth", &SpatRaster::setDepth)
		.method("setUnit", &SpatRaster::setUnit)
		.method("set_resolution", &SpatRaster::setResolution)
		.method("subset", &SpatRaster::subset)

		.method("cellFromXY", ( std::vector<double> (SpatRaster::*)(std::vector<double>,std::vector<double>) )( &SpatRaster::cellFromXY ))
		.method("vectCells", &SpatRaster::vectCells)
		.method("extCells", &SpatRaster::extCells)

		.method("wmean_rast", (SpatRaster (SpatRaster::*)(SpatRaster, bool, SpatOptions&))( &SpatRaster::weighted_mean ))
		.method("wmean_vect", (SpatRaster (SpatRaster::*)(std::vector<double>, bool, SpatOptions&))( &SpatRaster::weighted_mean ))

		.method("cellFromRowCol", ( std::vector<double> (SpatRaster::*)(std::vector<int_64>,std::vector<int_64>) )( &SpatRaster::cellFromRowCol ))
		.method("cellFromRowColCombine", ( std::vector<double> (SpatRaster::*)(std::vector<int_64>,std::vector<int_64>) )( &SpatRaster::cellFromRowColCombine ))
		.method("yFromRow", ( std::vector<double> (SpatRaster::*)(const std::vector<int_64>&) )( &SpatRaster::yFromRow ))
		.method("xFromCol", ( std::vector<double> (SpatRaster::*)(const std::vector<int_64>&) )( &SpatRaster::xFromCol ))
		.method("colFromX", ( std::vector<int_64> (SpatRaster::*)(const std::vector<double>&) )( &SpatRaster::colFromX ))
		.method("rowFromY", ( std::vector<int_64> (SpatRaster::*)(const std::vector<double>&) )( &SpatRaster::rowFromY ))
		.method("xyFromCell", ( std::vector< std::vector<double> > (SpatRaster::*)(std::vector<double>&) )( &SpatRaster::xyFromCell ))
		.method("rowColFromCell", ( std::vector< std::vector<int_64> > (SpatRaster::*)(std::vector<double>) )( &SpatRaster::rowColFromCell ))
		.method("readStart", &SpatRaster::readStart)
		.method("readStop", &SpatRaster::readStop)
		.method("readAll", &SpatRaster::readAll)
		.method("readValues", &SpatRaster::readValuesR)
		.method("getValues", &SpatRaster::getValues)
		.method("getBlockSizeR", &getBlockSizeR)
		.method("getBlockSizeWrite", &getBlockSizeWrite)
		.method("same", &sameObject)
#ifdef useRcpp
		.method("setValuesRcpp", &SpatRaster::setValuesRcpp)
#endif
		.method("setValues", &SpatRaster::setValues)
		.method("replaceCellValues", &SpatRaster::replaceCellValues)
		.method("replaceCellValuesLayer", &SpatRaster::replaceCellValuesLayer)
		.method("setRange", &SpatRaster::setRange)
		.method("writeStart", &SpatRaster::writeStart)
		.method("writeStop", &SpatRaster::writeStop)
		.method("writeValues", &SpatRaster::writeValues)
		.method("writeRaster", &SpatRaster::writeRaster)
		.method("canProcessInMemory", &SpatRaster::canProcessInMemory)
		.method("chunkSize", &SpatRaster::chunkSize)
//		.method("to_memory", &SpatRaster::to_memory, "to_memory")
		.method("update_meta", &SpatRaster::update_meta)


		.method("adjacentMat", &SpatRaster::adjacentMat)
		.method("adjacent", &SpatRaster::adjacent)
		.method("aggregate", &SpatRaster::aggregate)
		.method("align", &SpatRaster::align)
		.method("apply", &SpatRaster::apply)
		.method("rapply", &SpatRaster::rapply)
		.method("rappvals", &SpatRaster::rappvals)
		.method("roll", &SpatRaster::roll)
		.method("fill_range", &SpatRaster::fill_range)
		.method("arith_rast", ( SpatRaster (SpatRaster::*)(SpatRaster, std::string, SpatOptions&) )( &SpatRaster::arith ))
		.method("arith_numb", ( SpatRaster (SpatRaster::*)(std::vector<double>, std::string, bool, SpatOptions&) )( &SpatRaster::arith ))
		.method("rst_area", &SpatRaster::rst_area)
		.method("sum_area", &SpatRaster::sum_area)
		.method("area_by_value", &SpatRaster::area_by_value)

		.method("as_points", &SpatRaster::as_points)
		.method("as_points_value", &SpatRaster::as_points_value)
		.method("cells_notna", &SpatRaster::cells_notna)
		.method("cells_notna_novalues", &SpatRaster::cells_notna_novalues)
		.method("as_multipoints", &SpatRaster::as_multipoints)

		.method("as_lines", &SpatRaster::as_lines)
		.method("as_polygons", &SpatRaster::as_polygons)
		.method("polygonize", &SpatRaster::polygonize)

		.method("atan2", &SpatRaster::atan_2)

		.method("bilinearValues", &SpatRaster::bilinearValues)

		.method("patches", &SpatRaster::clumps)
		.method("boundaries", &SpatRaster::edges)
		.method("buffer", &SpatRaster::buffer)
		.method("gridDistance", &SpatRaster::gridDistance)
		.method("costDistance", &SpatRaster::costDistance)
		.method("rastDistance", &SpatRaster::distance)
		.method("rastDirection", &SpatRaster::direction)
		.method("make_tiles", &SpatRaster::make_tiles)
		.method("make_tiles_vect", &SpatRaster::make_tiles_vect)
		.method("ext_from_rc", &SpatRaster::ext_from_rc)

		.method("combineCats", &SpatRaster::combineCats)
		.method("droplevels", &SpatRaster::dropLevels)
		.method("vectDistanceDirect", &SpatRaster::distance_spatvector)
		.method("vectDistanceRasterize", &SpatRaster::distance_rasterize)
		.method("vectDirectionRasterize", &SpatRaster::direction_rasterize)
		.method("clamp", &SpatRaster::clamp)
		.method("replaceValues", &SpatRaster::replaceValues)
		.method("classify", ( SpatRaster (SpatRaster::*)(std::vector<double>, unsigned, unsigned, bool, bool, double, bool, bool, bool, SpatOptions&) )( &SpatRaster::reclassify))
		//.method("source_collapse", &SpatRaster::collapse, "collapse")
		.method("selRange", &SpatRaster::selRange)
		.method("separate", &SpatRaster::separate)
		.method("sort", &SpatRaster::sort)
		.method("intersect", &SpatRaster::intersect)

		.method("cover", &SpatRaster::cover)
		.method("crop", &SpatRaster::crop)
		.method("crop_mask", &SpatRaster::cropmask)
		.method("cum", &SpatRaster::cum)
		.method("disaggregate", &SpatRaster::disaggregate)
		.method("expand", &SpatRaster::extend)
		.method("extractCell", &SpatRaster::extractCell)
		.method("extractXY", &SpatRaster::extractXY)
//		.method("extractXYFlat", &SpatRaster::extractXYFlat, "extractXYflat")
		.method("extractVector", &SpatRaster::extractVector)
		.method("extractVectorFlat", &SpatRaster::extractVectorFlat)
		.method("flip", &SpatRaster::flip)
		.method("focal", &SpatRaster::focal)
		.method("focalValues", &SpatRaster::focal_values)
		.method("count", &SpatRaster::count)
		.method("freq", &SpatRaster::freq)
		.method("geometry", &SpatRaster::geometry)

		.method("get_aggregates", &SpatRaster::get_aggregates)
		.method("get_aggregate_dims", &SpatRaster::get_aggregate_dims2)
		.method("global", &SpatRaster::global)
		.method("global_weighted_mean", &SpatRaster::global_weighted_mean)

		.method("initf", ( SpatRaster (SpatRaster::*)(std::string, bool, SpatOptions&) )( &SpatRaster::init ), "init fun")
		.method("initv", ( SpatRaster (SpatRaster::*)(std::vector<double>, SpatOptions&) )( &SpatRaster::init ), "init value")
		.method("is_in", &SpatRaster::is_in)
		.method("is_in_cells", &SpatRaster::is_in_cells)
		.method("anynan", &SpatRaster::anynan)
		.method("nonan", &SpatRaster::nonan)
		.method("allnan", &SpatRaster::allnan)
		.method("isnan", &SpatRaster::isnan)
		.method("not_na", &SpatRaster::isnotnan)
		.method("isfinite", &SpatRaster::isfinite)
		.method("isinfinite", &SpatRaster::isinfinite)
		.method("logic_rast", ( SpatRaster (SpatRaster::*)(SpatRaster, std::string, SpatOptions&) )( &SpatRaster::logic ))
		.method("logic_numb", ( SpatRaster (SpatRaster::*)(bool, std::string, SpatOptions&) )( &SpatRaster::logic ))
		.method("mask_raster", ( SpatRaster (SpatRaster::*)(SpatRaster, bool, std::vector<double>, double, SpatOptions&) )( &SpatRaster::mask))
		.method("mask_vector", ( SpatRaster (SpatRaster::*)(SpatVector, bool, double, bool, SpatOptions&) )( &SpatRaster::mask))
		.method("math", &SpatRaster::math)
		.method("math2", &SpatRaster::math2)
		.method("modal", &SpatRaster::modal)
		.method("quantile", &SpatRaster::quantile)
		.method("rasterize", &SpatRaster::rasterize)
		.method("rasterizePoints", &SpatRaster::rasterizePoints)
		.method("rasterizeLyr", &SpatRaster::rasterizeLyr)
		.method("rasterizeGeom", &SpatRaster::rasterizeGeom)
		.method("rasterizeWindow", &SpatRaster::rasterizeWindow)
		.method("wincircle", &SpatRaster::win_circle)
		.method("winrect", &SpatRaster::win_rect)
		.method("rgb2col", &SpatRaster::rgb2col)
		.method("rgb2hsx", &SpatRaster::rgb2hsx)
		.method("hsx2rgb", &SpatRaster::hsx2rgb)

		.method("reverse", &SpatRaster::reverse)
		.method("rotate", &SpatRaster::rotate)
		//.method("sampleCells", &SpatRaster::sampleCells, "sampleCells")
		.method("sampleRegularRaster", &SpatRaster::sampleRegularRaster)
		.method("sampleRowColRaster", &SpatRaster::sampleRowColRaster)
		.method("sampleRegularValues", &SpatRaster::sampleRegularValues)
		.method("sampleRowColValues", &SpatRaster::sampleRowColValues)
		.method("sampleRandomRaster", &SpatRaster::sampleRandomRaster)
		.method("sampleRandomValues", &SpatRaster::sampleRandomValues)
		.method("scale", &SpatRaster::scale)
		.method("shift", &SpatRaster::shift)
		.method("terrain", &SpatRaster::terrain)
		.method("summary", &SpatRaster::summary)
		.method("summary_numb", &SpatRaster::summary_numb)
		.method("transpose", &SpatRaster::transpose)
		.method("trig", &SpatRaster::trig)
		.method("trim1", &SpatRaster::trim1)
		.method("trim", &SpatRaster::trim2)
		.method("unique", &SpatRaster::unique)
		.method("where", &SpatRaster::where)
		.method("sieve", &SpatRaster::sieveFilter)
		.method("view", &SpatRaster::viewshed)
		.method("proximity", &SpatRaster::proximity)

		.method("rectify", &SpatRaster::rectify)
		.method("stretch", &SpatRaster::stretch)
		.method("warp", &SpatRaster::warper)
		.method("resample", &SpatRaster::resample)
		.method("zonal", &SpatRaster::zonal)
		.method("zonal_old", &SpatRaster::zonal_old)
		.method("is_true", &SpatRaster::is_true)
		.method("is_false", &SpatRaster::is_false)
	;

    class_<SpatRasterCollection>("SpatRasterCollection")
		.constructor()
	    .constructor<std::string, std::vector<int>, bool>()

		.property("names", &SpatRasterCollection::get_names, &SpatRasterCollection::set_names)

		.method("deepcopy", &SpatRasterCollection::deepCopy)
		.method("dims", &SpatRasterCollection::dims)
		.method("extent", &SpatRasterCollection::getExtent)

		.method("has_error", &SpatRasterCollection::has_error)
		.method("has_warning", &SpatRasterCollection::has_warning)
		.method("getError", &SpatRasterCollection::getError)
		.method("getWarnings", &SpatRasterCollection::getWarnings)
		//.field("messages", &SpatRasterCollection::msg, "messages")
		.field_readonly("x", &SpatRasterCollection::ds)
		.method("length", &SpatRasterCollection::size)
		.method("resize", &SpatRasterCollection::resize)
		.method("erase", &SpatRasterCollection::erase)
		.method("add", &SpatRasterCollection::push_back)
		.method("merge", &SpatRasterCollection::merge)
		.method("mosaic", &SpatRasterCollection::mosaic)
		.method("morph", &SpatRasterCollection::morph)
	;

    class_<SpatRasterStack>("SpatRasterStack")
		.constructor()
	    .constructor<std::string, std::vector<int>, bool>()
	    .constructor<SpatRaster, std::string, std::string, std::string>()
		.method("deepcopy", &SpatRasterStack::deepCopy)

		.method("has_error", &SpatRasterStack::has_error)
		.method("has_warning", &SpatRasterStack::has_warning)
		.method("getError", &SpatRasterStack::getError)
		.method("getWarnings", &SpatRasterStack::getWarnings)

		.method("readStart", &SpatRasterStack::readStart)
		.method("readStop", &SpatRasterStack::readStop)
		.method("nsds", &SpatRasterStack::nsds)
		.method("ncol", &SpatRasterStack::ncol)
		.method("nrow", &SpatRasterStack::nrow)
		.method("nlyr", &SpatRasterStack::nlyr)
		.method("res", &SpatRasterStack::resolution)
		.method("ext", &SpatRasterStack::getExtent)
		.method("filenames", &SpatRasterStack::filenames)
		
		.method("get_crs", &SpatRasterStack::getSRS)
		.property("names", &SpatRasterStack::get_names, &SpatRasterStack::set_names)
		.property("long_names", &SpatRasterStack::get_longnames, &SpatRasterStack::set_longnames)
		.property("units", &SpatRasterStack::get_units, &SpatRasterStack::set_units)
		.method("add", &SpatRasterStack::push_back)
		.method("resize", &SpatRasterStack::resize)
		.method("summary", &SpatRasterStack::summary)
		.method("summary_numb", &SpatRasterStack::summary_numb)
		.method("getsds", &SpatRasterStack::getsds)
		.method("replace", &SpatRasterStack::replace)
		.method("subset", &SpatRasterStack::subset)
		.method("collapse", &SpatRasterStack::collapse )
		.method("extractCell", &SpatRasterStack::extractCell)
		.method("extractVector", &SpatRasterStack::extractVector)
		.method("crop", &SpatRasterStack::crop)
	;
}


/*
SpatRaster SQRT() {
	SpatRaster r = *this;
	std::transform(r.values.begin(), r.values.end(), r.values.begin(), (double(*)(double)) sqrt);
	return r;
}

SpatRaster SQRTfree(SpatRaster* g) {
	SpatRaster r = *g;
	std::transform(r.values.begin(), r.values.end(), r.values.begin(), (double(*)(double)) sqrt);
	return r;
}
*/

