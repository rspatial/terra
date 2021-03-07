#include "spatRaster.h"

#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"

#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"


bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {

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

	bool warngroup = false;
	std::vector<std::string> gnames;
	if (sub == "") {
		char** papszOptions = NULL;
		gnames = poRootGroup->GetMDArrayNames(papszOptions);
		CSLDestroy(papszOptions);
		sub = gnames[0];
		if (gnames.size() > 1) warngroup = true;
	}

    auto poVar = poRootGroup->OpenMDArray(sub.c_str());
    if( !poVar )   {
		setError("cannot find: " + sub);
		return false;
    }

	SpatRasterSource s;

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
	std::string msg;
	if (!s.srs.set({wkt}, msg)) {
		addWarning(msg);
	}


	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames;
    for( const auto &poDim: poVar->GetDimensions() ) {
        dimcount.push_back(static_cast<size_t>(poDim->GetSize()));
        dimnames.push_back(static_cast<std::string>(poDim->GetName()));
    }
	s.m_ndims = dimcount.size();

	if (warngroup) {
		std::string gn = "";
		for (size_t i=0; i<gnames.size(); i++) {
			if (!is_in_vector(gnames[i], dimnames) &&  (gnames[i] != sub)) { 
				gn += gnames[i] + ", ";
			}
		}
		addWarning("using: " + sub + ". Other groups are: \n" + gn);
	}
	s.source_name = sub;
	s.source_name_long = poVar->GetAttribute("long_name")->ReadAsString();

	s.m_hasNA = false;
	double NAval = poVar->GetNoDataValueAsDouble(&s.m_hasNA);
	if (s.m_hasNA) {
		s.m_missing_value = NAval;
		Rcpp::Rcout << "NAval: " << NAval << std::endl;
	}

	if (xyz[0] < s.m_ndims) {
		s.ncol = dimcount[xyz[0]];
		s.m_dimnames.push_back(dimnames[xyz[0]]);
	} else {
		setError("the first dimension is not valid");
		return false;		
	}
	if (xyz[1] < s.m_ndims) {
		s.nrow = dimcount[xyz[1]];
		s.m_dimnames.push_back(dimnames[xyz[1]]);
	} else {
		setError("the second dimension is not valid");
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
	
	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);

	std::vector<GUInt64> offset(source[src].m_ndims, 0);
	offst[dims[0]] = row;
	offst[dims[1]] = col;

//	std::vector<size_t> count = source[src].m_counts;
	std::vector<size_t> count(source[src].m_ndims, 1);
	count[dims[0]] = nrows;
	count[dims[1]] = ncols;
	count[dims[2]] = nlyr();

	size_t n=1;
	//count = {3600, 1, 1, 7200, 1};
	for (size_t i=0; i<count.size(); i++) {
		//count[i] = std::min(count[i], source[src].m_counts[i]);
		Rcpp::Rcout << offset[i] << "-" << count[i] << ", ";
		n *= count[i];
	}
	Rcpp::Rcout << std::endl;
	
	out.resize(n, -99);
		
    GDALMDArrayRead(source[src].gdalmdarray,
                    &offset[0],
                    &count[0],
                    NULL, /* step: defaults to 1,1,1 */
                    NULL, /* stride: default to row-major convention */
                    hDT,
                    &out[0],
                    NULL, /* array start. Omitted */
                    0 /* array size in bytes. Omitted */);
    GDALExtendedDataTypeRelease(hDT);
 	
	std::replace (out.begin(), out.end(), source[src].m_missing_value, (double)NAN);
	
	return true;
}

  
#else  


bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {
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

#endif

