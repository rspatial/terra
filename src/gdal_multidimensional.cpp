
#include "spatRaster.h"

/*
#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"

#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"


bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> sub, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<size_t> xyz) {

Rcpp::Rcout << "in" << std::endl;

	if (xyz.size() != 3) {
		setError("you must supply three dimension indices");
        return false;
	}


    auto poDataset = std::unique_ptr<GDALDataset>(
        GDALDataset::Open(fname.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !poDataset ) {
		setError("not a good dataset");
        return false;
    }

	auto poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no roots");
		return false;
    }


	std::vector<std::string> gnames;
	std::string subdsname = "";
	char** papszOptions = NULL;
	gnames = poRootGroup->GetMDArrayNames(papszOptions);
	CSLDestroy(papszOptions);

	if (gnames.size() == 0) {
		setError("no subdatsets detected");
		return false;
	}
	Rcpp::Rcout << "available: ";
	for (size_t i=0; i<gnames.size(); i++) {
		Rcpp::Rcout << gnames[i] << " ";
	}
	Rcpp::Rcout << std::endl;


	if (subname.size() > 0) {
		subdsname = subname[0];
		if (std::find(gnames.begin(), gnames.end(), subdsname) == gnames.end()) {
			setError("subdatset name not found");
			return false;
		}
	} else if (sub[0] >= 0) {
		if (sub[0] >= (int)gnames.size()) {
			setError("subdatset is out or range");
			return false;
		} else {
			subdsname = gnames[sub[0]];
		}
	} else {
		subdsname = gnames[0];
		if (gnames.size() > 1)  {
			std::string gn = "";
			for (size_t i=1; i<gnames.size(); i++) {
				gn += gnames[i] + ", ";
			}
			addWarning("using: " + subdsname + ". Other groups are: \n" + gn);
		}
	}

	Rcpp::Rcout << "subdsname: " << subdsname << std::endl;

    auto poVar = poRootGroup->OpenMDArray(subdsname.c_str());
    if( !poVar )   {
		setError("cannot find: " + subdsname);
		return false;
    }


	std::string wkt = "";
	std::shared_ptr<OGRSpatialReference> srs = poVar->GetSpatialRef();
	if (srs != NULL) {
		char *cp;
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = srs->exportToWkt(&cp, options);
		if (err == OGRERR_NONE) {
			wkt = std::string(cp);
		}
		CPLFree(cp);
	}

	SpatRasterSource s;

	std::string msg;
	if (!s.srs.set({wkt}, msg)) {
		addWarning(msg);
	}

	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames;
	std::vector<double> dim_start, dim_end;

	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);
    for ( const auto &poDim: poVar->GetDimensions() ) {
        dimcount.push_back(static_cast<size_t>(poDim->GetSize()));
        dimnames.push_back(static_cast<std::string>(poDim->GetName()));

		std::vector<size_t> count = {dimcount[dimcount.size()-1]};
		std::vector<double> vals(count[0]);

		const auto indvar = poDim->GetIndexingVariable();
		indvar->Read(
			std::vector<GUInt64>{0}.data(),
            count.data(), nullptr, nullptr,
            GDALExtendedDataType::Create(GDT_Float64),
            &vals[0]);

		// to do: check for equal spacing if x or y dim
		dim_start.push_back(vals[0]);
        dim_end.push_back(vals[vals.size()-1]);
		Rcpp::Rcout << vals[0] << " - " << vals[vals.size()-1] << std::endl;

    }
    GDALExtendedDataTypeRelease(hDT);

	s.m_ndims = dimcount.size();
	s.source_name = subdsname;
	s.source_name_long = poVar->GetAttribute("long_name")->ReadAsString();

	s.m_hasNA = false;
	double NAval = poVar->GetNoDataValueAsDouble(&s.m_hasNA);
	if (s.m_hasNA) {
		s.m_missing_value = NAval;
	}

	SpatExtent e;
	if (xyz[0] < s.m_ndims) {
		s.nrow = dimcount[xyz[0]];
		s.m_dimnames.push_back(dimnames[xyz[0]]);
		double res = (dim_start[xyz[0]] - dim_end[xyz[0]]) / (s.nrow-1);
		e.ymax = dim_start[xyz[0]] + 0.5 * res;
		e.ymin = dim_end[xyz[0]] - 0.5 * res;
	} else {
		setError("the second dimension is not valid");
		return false;
	}
	if (xyz[1] < s.m_ndims) {
		s.ncol = dimcount[xyz[1]];
		s.m_dimnames.push_back(dimnames[xyz[1]]);
		double res = (dim_end[xyz[1]] - dim_start[xyz[1]]) / (s.ncol-1);
		e.xmin = dim_start[xyz[1]] - 0.5 * res;
		e.xmax = dim_end[xyz[1]] + 0.5 * res;
	} else {
		setError("the first dimension is not valid");
		return false;
	}
	if (s.m_ndims > 2) {
		if (xyz[2] < s.m_ndims) {
			s.nlyr = dimcount[xyz[2]];
			s.m_dimnames.push_back(dimnames[xyz[2]]);
		} else {
			setError("the third dimension is not valid");
			return false;
		}
	}
	s.m_dims = xyz;
	s.extent = e;
	if (s.m_ndims > 3) {
		for (size_t i=0; i<s.m_ndims; i++) {
			bool found = false;
			for (size_t j=0; j<3; j++) {
				if (i == xyz[j]) found = true;
			}
			if (!found) {
				s.m_dims.push_back(i);
				s.m_dimnames.push_back(dimnames[i]);
			}
		}
	}

	s.nlyrfile = s.nlyr;
	s.resize(s.nlyr);
	s.layers.resize(s.nlyr);
    std::iota(s.layers.begin(), s.layers.end(), 0);
	s.flipped = false;
	s.rotated = false;
	s.memory = false;
	s.filename = fname;
	s.hasValues = true;
	s.unit = std::vector<std::string>(s.nlyr, poVar->GetUnit());
	s.multidim = true;

// layer names
	//std::vector<std::string> nms(s.nlyr, "");
	//s.names = nms;
// time
// extent
	s.m_counts = dimcount;
	setSource(s);
	for (size_t i=0; i<s.m_ndims; i++){
		Rcpp::Rcout << s.m_dims[i] << " " << s.m_dimnames[i] << " " << s.m_counts[s.m_dims[i]] << std::endl;
	}
	return true;
}



bool SpatRaster::readStartMulti(unsigned src) {

    GDALDatasetH hDS = GDALOpenEx( source[src].filename.c_str(), GDAL_OF_MULTIDIM_RASTER, NULL, NULL, NULL);
    if (!hDS) {
		setError("not a good dataset");
        return false;
    }
    GDALGroupH hGroup = GDALDatasetGetRootGroup(hDS);
    GDALReleaseDataset(hDS);
    if (!hGroup) {
		setError("not a good root group");
		return false;
    }

	GDALMDArrayH hVar = GDALGroupOpenMDArray(hGroup, source[src].source_name.c_str(), NULL);
    GDALGroupRelease(hGroup);
    if (!hVar) {
		setError("not a good array");
		return false;
    }

	source[src].gdalmdarray = hVar;
	return true;
}


bool SpatRaster::readStopMulti(unsigned src) {
	GDALMDArrayRelease(source[src].gdalmdarray);
	source[src].open_read = false;
	return true;
}


bool SpatRaster::readValuesMulti(std::vector<double> &out, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {

	Rcpp::Rcout << "reading" << std::endl;

	std::vector<GUInt64> offset(source[src].m_ndims, 0);
	std::vector<size_t> dims = source[src].m_dims;

	offset[source[src].m_dims[0]] = col;
	offset[source[src].m_dims[1]] = row;
	offset[source[src].m_dims[2]] = 0;

//	std::vector<size_t> count = source[src].m_counts;
	std::vector<size_t> count(source[src].m_ndims, 1);
	count[source[src].m_dims[0]] = ncols;
	count[source[src].m_dims[1]] = nrows;
	count[source[src].m_dims[2]] = nlyr();

	size_t n=1;
	for (size_t i=0; i<count.size(); i++) {
		n *= count[i];
	}

	//count = {3600, 1, 1, 7200, 1};

	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);

	std::vector<double> temp;
	temp.resize(n);
	GDALMDArrayRead(source[src].gdalmdarray,
						&offset[0],
						&count[0],
						NULL, // step: defaults to 1,1,1
						NULL, // stride: default to row-major convention
						hDT,
						&temp[0],
						NULL, // array start. Omitted
						0 // array size in bytes. Omitted
						);
    GDALExtendedDataTypeRelease(hDT);


//tbd: row order should be reversed

//	size_t nc = nrows * ncols;
//	size_t nl = nlyr();
//	out.resize(0);
//	out.reserve(n);
//	for (size_t i=0; i<nl; i++) {
//		for (size_t j=0; j<nc; j++) {
//			out.push_back( temp[nl*j + i] );
//		}
//	}


    out = std::move(temp);
	if (source[src].m_hasNA) {
		std::replace (out.begin(), out.end(), source[src].m_missing_value, (double)NAN);
	}

	return true;
}


#else

*/

bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> sub, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<size_t> xyz) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}

bool SpatRaster::readStartMulti(unsigned src) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}
bool SpatRaster::readStopMulti(unsigned src) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}


bool SpatRaster::readValuesMulti(std::vector<double> &out, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}


//#endif

