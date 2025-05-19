
#include "spatRaster.h"


#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 1

#include "proj.h"

#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"
#include "vecmath.h"



bool parse_ncdf_time(SpatRasterSource &s, const std::string unit, const std::string calendar, std::vector<double> raw, std::string &msg) {

	std::vector<int_64> out;
	std::string origin = unit;
	bool years = false;
	bool yearsbp = false;
	bool yearmonths = false;
	bool months = false;
	bool days = false;
	bool hours = false;
	bool minutes = false;
	bool seconds = false;
	bool foundorigin = false;
	std::string step;
	
	lowercase(origin);
	if ((origin.find("seconds")) != std::string::npos) {
		seconds = true;
	} else if ((origin.find("minutes")) != std::string::npos) {
		minutes = true;
	} else if ((origin.find("hours")) != std::string::npos) {
		hours = true;
	} else if ((origin.find("days")) != std::string::npos) {
		days = true;
	} else if ((origin.find("months since")) != std::string::npos) {
		yearmonths = true;
		foundorigin = true;
	} else if ((origin.find("months")) != std::string::npos) {
		months = true;
		foundorigin = true;
	} else if ((origin.find("years before present")) != std::string::npos) {
		yearsbp = true;
		foundorigin = true;
	} else if ((origin.find("years")) != std::string::npos) {
		years = true;
	}

	if (!foundorigin) {
		size_t pos;		
		if ((pos = origin.find("from")) != std::string::npos) {
			origin.erase(0, pos + 5);
			foundorigin = true;
		} else if ((pos = origin.find("since")) != std::string::npos) {
			origin.erase(0, pos + 6);
			foundorigin = true;
		}
	}

	SpatTime_t offset = 0;
	if (foundorigin) {
		step = "seconds";
		out.reserve(raw.size());
		std::string cal = "366";
		if (calendar == "360_day" || calendar == "360 day") {
			cal = "360";
		} else if (calendar == "noleap" || calendar == "365_day" || calendar == "365 day") {
			cal = "365";		
		} else if (calendar =="gregorian" || calendar =="proleptic_gregorian" || calendar=="standard" || calendar == "julian") {
			cal = "366";
		} else if (!(months || years || yearmonths || yearsbp)) {
			//cal = "366";
			msg = "unknown calendar (assuming standard): " + calendar;			
		}

		// this shortcut means that 360/noleap calendars loose only have dates, no time
		// to be refined
		if ((hours || minutes || seconds) && (cal == "360")) {
			int div = 24;
			double add = 0;
			std::vector<int> ymd = getymd(origin);
			if (hours) {
				hours = false;
				add = ymd[3] + ymd[4] / 60 + ymd[5] / 3600; 
			} else if (minutes) {
				div = 1440; // 24 * 60
				add = ymd[3] * 60 + ymd[4] + ymd[5] / 60; 
				minutes = false;
			} else if (seconds) {
				div = 86400; // 24 * 3600
				add = ymd[3] * 3600 + ymd[4] * 60 + ymd[5]; 
				seconds = false;
			}
			for (size_t i=0; i<raw.size(); i++) {
				raw[i] = (raw[i]+add) / div; 
			}
			days = true;
		} 
	
		if (days) {
			step = "days";
			std::vector<int> ymd = getymd(origin);
			if (cal == "365") {
				for (size_t i=0; i<raw.size(); i++) out.push_back(
					get_time_noleap(ymd[0], ymd[1], ymd[2], 0, 0, 0, raw[i], "days"));
			} else if (cal == "360") {
				for (size_t i=0; i<raw.size(); i++) out.push_back(
					time_from_day_360(ymd[0], ymd[1], ymd[2], raw[i]));
			} else {
				for (size_t i=0; i<raw.size(); i++) out.push_back(
					time_from_day(ymd[0], ymd[1], ymd[2], raw[i]));
			}
		} else if (hours) {
			if (cal == "365") {
				std::vector<int> ymd = getymd(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(
						get_time_noleap(ymd[0], ymd[1], ymd[2], ymd[3], 0, 0, raw[i], "hours")
					);
			} else {
				offset = get_time_string(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(raw[i]*3600+offset);
			}
		} else if (minutes) {
			if (cal == "365") {
				std::vector<int> ymd = getymd(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(
						get_time_noleap(ymd[0], ymd[1], ymd[2], ymd[3], ymd[4], 0, raw[i], "minutes")
					);
			} else {
				offset = get_time_string(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(60*raw[i]+offset);
			}
		} else if (seconds) {
			if (cal == "365") {
				std::vector<int> ymd = getymd(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(
						get_time_noleap(ymd[0], ymd[1], ymd[2], ymd[3], ymd[4], 0, raw[i], "minutes")
					);
			} else {
				offset = get_time_string(origin);
				for (size_t i=0; i<raw.size(); i++) out.push_back(raw[i]+offset);
			}
		} else if (years) {
			step = "years";
			int syear = getyear(origin);
			for (size_t i=0; i<raw.size(); i++) out.push_back(get_time(syear+raw[i], 6, 30, 0, 0, 0));
		} else if (yearsbp) {
			step = "years";
			int syear = 1950;
			for (size_t i=0; i<raw.size(); i++) out.push_back(get_time(syear+raw[i], 6, 30, 0, 0, 0));
		} else if (yearmonths) {
			step = "yearmonths";
			int syear = getyear(origin);
			for (size_t i=0; i<raw.size(); i++) {
				long year = std::floor(raw[i] / 12.0);
				int month = raw[i] - (12 * year) + 1;
				out.push_back(get_time(syear+year, month, 15, 0, 0, 0));
			}
		} else if (months) {
			step = "months";
			// check for 0..11 range
			int zero = vmin(raw, true) == 0.0;
			for (size_t i=0; i<raw.size(); i++) {
				unsigned m = std::ceil(raw[i] + zero);
				out.push_back(get_time(1970, m, 15, 0, 0, 0));
			}
		} else {
			step = "raw";
			for (size_t i=0; i<raw.size(); i++) out.push_back(raw[i]);
		}
	}

	s.time = out;
	s.timestep = step;
	s.hasTime = true;
	return true;
}





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

	bool verbose = false;

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

	
	if (!dimfo(poRootGroup, names, dim_names, dim_size, msg)) {
		setError(msg);
		return false;
	} else {
		for (size_t i=0; i<names.size(); i++) {
			size_t ni = dim_size[i].size();
			if (ni > 1) {
				ar_names.push_back(names[i]);
				if (verbose) {
					Rcpp::Rcout << names[i] << ": ";	
					for (size_t j=0; j<ni; j++) {
						Rcpp::Rcout << dim_names[i][j] << " (" << dim_size[i][j] << ") ";	
					}
					Rcpp::Rcout << std::endl;
				}
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
	} 
	msg = "";
	if (!s.srs.set({wkt}, msg)) {
		addWarning(msg);
	}


//	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);
//  GDALExtendedDataTypeRelease(hDT);

// dimensions 
	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames, dimunits;
//	std::vector<double> dim_start, dim_end;

	std::vector<std::vector<double>> dimvals;
	dimvals.reserve(dim_names[sub[0]].size());

	size_t i=0;
    for ( const auto &poDim: poVar->GetDimensions() ) {
		size_t n = poDim->GetSize();
        dimcount.push_back(static_cast<size_t>(n));
        dimnames.push_back(static_cast<std::string>(poDim->GetName()));
			
		std::vector<GUInt64> start(1, 0);
		std::vector<size_t> count = {n};
		dimvals.push_back(std::vector<double>(n));

		const auto indvar = poDim->GetIndexingVariable();
        dimunits.push_back(static_cast<std::string>(indvar->GetUnit()));

		indvar->Read(start.data(), count.data(), nullptr, nullptr,  GDALExtendedDataType::Create(GDT_Float64), &dimvals[i][0]);
		i++;
    }

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
	s.nrow = dimcount[ix];
	s.nlyr = 1;

	// to do: check for equal spacing if x or y dim

	double ystart = dimvals[iy][dimvals[iy].size()-1];
	double yend = dimvals[iy][0];
	double res = (yend - ystart) / (s.ncol-1);
	e.ymax = yend + 0.5 * res;
	e.ymin = ystart - 0.5 * res;
	s.m_dimnames.push_back(dimnames[iy]);

	double xstart = dimvals[ix][0];
	double xend = dimvals[ix][dimvals[ix].size()-1];
	res = (xend - xstart) / (s.nrow-1);
	e.xmin = dimvals[ix][0] - 0.5 * res;
	e.xmax = dimvals[ix][dimvals[ix].size()-1] + 0.5 * res;

	s.m_dimnames.push_back(dimnames[ix]);

	if (s.m_ndims > 2) {
		s.nlyr = dimcount[it];
		s.m_dimnames.push_back(dimnames[it]);
		std::string msg;
		parse_ncdf_time(s, dimunits[it], "standard", dimvals[it], msg);
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

	std::vector<std::string> nms(s.nlyr, subdsname + "_");
	for (size_t i=0; i<=nms.size(); i++) {
		nms[i] += std::to_string(i + 1);
	}	
	s.names = nms;

// time

	s.m_counts = dimcount;
	setSource(s);
	if (verbose) {
		for (size_t i=0; i<s.m_ndims; i++){
			Rcpp::Rcout << s.m_dims[i] << " " << s.m_dimnames[i] << " " << s.m_counts[i] << std::endl;
		}
	}
	return true;
}



bool SpatRaster::readStartMulti(size_t src) {

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
	GDALMDArrayRelease(source[src].gdalmdarray);
	source[src].open_read = false;
	return true;
}


bool SpatRaster::readValuesMulti(std::vector<double> &out, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {

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
	GDALMDArrayRead(source[src].gdalmdarray, &offset[0], &count[0],
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

