#include "spatRaster.h"

#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"

#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"


bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {

	size_t ndims = xyz.size();
	if (ndims != 3) {
		setError("need three dimension indices");
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

	s.multidim = true;

	std::vector<size_t> dims;
	std::vector<std::string> dimnames;
    for( const auto poDim: poVar->GetDimensions() ) {
        dims.push_back(static_cast<size_t>(poDim->GetSize()));
        dimnames.push_back(static_cast<std::string>(poDim->GetName()));
    }
	s.m_ndims = dims.size();

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

	if (xyz[0] < s.m_ndims) {
		s.ncol = dims[xyz[0]];
		s.m_dimnames.push_back(dimnames[xyz[0]]);
	} else {
		setError("the first dimension is not valid");
		return false;		
	}
	if (xyz[1] < s.m_ndims) {
		s.nrow = dims[xyz[1]];
		s.m_dimnames.push_back(dimnames[xyz[1]]);
	} else {
		setError("the second dimension is not valid");
		return false;		
	}
	if (s.m_ndims > 2) {
		if (xyz[2] < s.m_ndims) {
			s.nlyr = dims[xyz[2]];
			s.m_dimnames.push_back(dimnames[xyz[2]]);
		} else {
			setError("the third dimension is not valid");
			return false;		
		}
	}
	s.m_dims = xyz;
		
	if (dims.size() > 3) {
		for (size_t i=0; i<ndims; i++) {
			bool found = false;
			for (size_t j=0; j<3; j++) {
				if (i == xyz[j]) found = true;
			}
			if (!found) {
				s.m_dims.push_back(dims[i]);
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

// layer names
// time 
// extent

	setSource(s);
	for (size_t i=0; i<s.m_dims.size(); i++){
		Rcpp::Rcout << s.m_dims[i] << " " << s.m_dimnames[i] << " " << dims[s.m_dims[i]] << std::endl;
	}
	return true;
}


/*
read
    std::vector<double> values(nValues);
    poVar->Read(std::vector<GUInt64>{0,0,0}.data(),
                anCount.data(),
                nullptr, // step: defaults to 1,1,1 
                nullptr, // stride: default to row-major convention
                GDALExtendedDataType::Create(GDT_Float64),
                &values[0]);

			
	return values;
*/

  
#else  


bool SpatRaster::constructFromFileMulti(std::string fname, std::string sub, std::vector<size_t> xyz) {
	setError("multidim is not supported by GDAL < 3.1");
	return false;
}


#endif

