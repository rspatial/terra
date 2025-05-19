
#include "spatRaster.h"


#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"

#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"



bool dimfo(std::shared_ptr<GDALGroup> poRootGroup, std::vector<std::string> &ar_names, std::vector<std::vector<std::string>> &dimnames, std::vector<std::vector<size_t>> &dimsize, std::string &msg) {

	msg = "";
	char** papszOptions = NULL;
	ar_names = poRootGroup->GetMDArrayNames(papszOptions);
	CSLDestroy(papszOptions);

	size_t n = ar_names.size();
	if (n == 0) {
		msg = "no arrays detected";
		return false;
	}
	dimnames.resize(n);
	dimsize.resize(n);
	
	for (size_t i=0; i<ar_names.size(); i++) {
		auto poVar = poRootGroup->OpenMDArray(ar_names[i].c_str());
		if( !poVar )   {
			msg = ("cannot open: " + ar_names[i]);
			return false;
		}
		for ( const auto &poDim: poVar->GetDimensions() ) {
			dimnames[i].push_back(static_cast<std::string>(poDim->GetName()));
			dimsize[i].push_back(static_cast<size_t>(poDim->GetSize()));
		}
	}
	return true;
}



bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> sub, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<size_t> xyz) {

	SpatRasterSource s;

    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(fname.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !poDataset ) {
		setError("not a good dataset");
        return false;
    }

	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no roots");
		return false;
    }

	std::vector<std::string> names, ar_names;
	std::vector<std::vector<std::string>> dim_names;
	std::vector<std::vector<size_t>> dim_size;
	std::string msg;

	Rcpp::Rcout << "--- info ---" << std::endl;
	
	if (!dimfo(poRootGroup, names, dim_names, dim_size, msg)) {
		setError(msg);
		return false;
	} else {
		for (size_t i=0; i<names.size(); i++) {
			size_t ni = dim_size[i].size();
			if (ni > 1) {
				ar_names.push_back(names[i]);
				Rcpp::Rcout << names[i] << ": ";	
				for (size_t j=0; j<ni; j++) {
					Rcpp::Rcout << dim_names[i][j] << " (" << dim_size[i][j] << ") ";	
				}
				Rcpp::Rcout << std::endl;
			}
		}
	}

//	if (xyz.size() != 3) {
//		setError("you must supply three dimension indices");
//       return false;
//	}

	std::string subdsname = "";
	if (subname.size() > 0) {
		subdsname = subname[0];
		int w = where_in_vector(subdsname, ar_names, false);
		if (w < 0) {
			setError("array " + subdsname + " not found. Should be one of:\n  " + concatenate(ar_names, ", "));
			return false;
		} else {
			sub = {w};
		} 
	} else if (sub[0] >= 0) {
		if (sub[0] >= (int)ar_names.size()) {
			setError("array number is out or range");
			return false;
		} else {
			subdsname = ar_names[sub[0]];
		}
	} else {
		sub = {0};
		subdsname = ar_names[0];
		if (ar_names.size() > 1)  {
			addWarning("using array: " + subdsname + ". Other groups are: \n" + concatenate(ar_names, ", "));
		}
	}

    auto poVar = poRootGroup->OpenMDArray(subdsname.c_str());
    if( !poVar )   {
		setError("cannot find: " + subdsname);
		return false;
    }


	std::string wkt = "";
	auto srs = poVar->GetSpatialRef();

	if (srs != NULL) {
		char *cp;
		const char *options[3] = { "MULTILINE=YES", "FORMAT=WKT2", NULL };
		OGRErr err = srs->exportToWkt(&cp, options);
		if (err == OGRERR_NONE) {
			wkt = std::string(cp);
		}
		CPLFree(cp);
	} else {
		Rcpp::Rcout << "wkt is null" <<std::endl;		
	}


	msg = "";
	if (!s.srs.set({wkt}, msg)) {
		addWarning(msg);
	}


//	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);
//  GDALExtendedDataTypeRelease(hDT);

	Rcpp::Rcout << "--- dimensions ---" << std::endl;

// dimensions 
	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames;
	std::vector<double> dim_start, dim_end;


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
	Rcpp::Rcout << "--- end dimensions ---" << std::endl;


	s.nlyr = 1;

	s.m_ndims = dimcount.size();
	s.source_name = subdsname;
	auto lname = poVar->GetAttribute("long_name");
	if (lname) s.source_name_long = lname->ReadAsString();
	
	s.m_hasNA = false;
	double NAval = poVar->GetNoDataValueAsDouble(&s.m_hasNA);
	if (s.m_hasNA) {
		s.m_missing_value = NAval;
	}

	if (xyz.size() < 2) {
		xyz = {3,2,1,0};
	}

	xyz.resize(dimcount.size());
	if (xyz.size() < 2) {
		setError("insufficient dimensions");
		return false;
	}
	

	int ix = dimcount.size()-1;
	int iy = ix - 1;
	int it = ix - 2;
//	int iz = ix - 3;

	
	SpatExtent e;
	
 	s.ncol = dimcount[iy];
	s.m_dimnames.push_back(dimnames[iy]);
	double res = (dim_start[iy] - dim_end[iy]) / (s.ncol-1);
	e.ymax = dim_start[iy] + 0.5 * res;
	e.ymin = dim_end[iy] - 0.5 * res;
	s.m_dimnames.push_back(dimnames[iy]);

	s.nrow = dimcount[ix];
	s.m_dimnames.push_back(dimnames[ix]);
	res = (dim_end[ix] - dim_start[ix]) / (s.nrow-1);
	e.xmin = dim_start[ix] - 0.5 * res;
	e.xmax = dim_end[ix] + 0.5 * res;
	s.m_dimnames.push_back(dimnames[ix]);

	if (s.m_ndims > 2) {
		s.nlyr = dimcount[it];
		s.m_dimnames.push_back(dimnames[it]);
	}
	s.m_dims = xyz;
	s.extent = e;

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



bool SpatRaster::readStartMulti(size_t src) {

Rcpp::Rcout <<  "readStartMulti\n";

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


bool SpatRaster::readStopMulti(size_t src) {

Rcpp::Rcout <<  "readStopMulti\n";

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


bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> sub, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<size_t> xyz) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}

bool SpatRaster::readStartMulti(size_t src) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}
bool SpatRaster::readStopMulti(size_t src) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}


bool SpatRaster::readValuesMulti(std::vector<double> &out, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}


#endif

