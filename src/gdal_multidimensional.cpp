
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


inline bool ncdf_keep(std::string const &s) {
	std::vector<std::string> end = {"_bnds", "_bounds", "lat", "lon", "longitude", "latitude"};
	for (size_t i=0; i<end.size(); i++) {
		if (s.length() >= end[i].length()) {
			if (s.compare(s.length() - end[i].length(), s.length(), end[i]) == 0) {
				return false;
			}
		}
	}
	if (s == "x" || s == "y" || s == "northing" || s == "easting") {
		return false;
	}
	return true;
}


bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<int> sub, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> xyz) {

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
			if ((ni > 1) && (ncdf_keep(names[i]))) {
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

// dimensions 
	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames, dimunits;
	std::vector<std::vector<double>> dimvals;
	dimvals.reserve(dim_names[sub[0]].size());

	std::string calendar = "";	
	std::vector<std::shared_ptr<GDALDimension>> dimData = poVar->GetDimensions();
	size_t ndim = dimData.size();
    for (size_t i=0; i<ndim; i++) {
		size_t n = dimData[i]->GetSize();
        dimcount.push_back(static_cast<size_t>(n));
        std::string name = static_cast<std::string>(dimData[i]->GetName());
		dimnames.push_back(name);
		std::vector<GUInt64> start(1, 0);
		std::vector<size_t> count = {n};
		dimvals.push_back(std::vector<double>(n));

		if (verbose) Rcpp::Rcout << name << std::endl;

		const auto indvar = dimData[i]->GetIndexingVariable();
        dimunits.push_back(static_cast<std::string>(indvar->GetUnit()));
		auto pcal = indvar->GetAttribute("calendar");
		if (pcal) calendar = pcal->ReadAsString();
		
	
		indvar->Read(start.data(), count.data(), nullptr, nullptr, GDALExtendedDataType::Create(GDT_Float64), &dimvals[i][0]);
		if ((i >= (ndim-2)) && (dimvals[i].size() > 2)) {
			double res = dimvals[i][1] - dimvals[i][0];
			if (!indvar->IsRegularlySpaced(dimvals[i][0], res)) {
				setError(name + " is not regularly spaced");
				return false;
			}
		}
    }

	s.m_ndims = dimcount.size();
	s.source_name = subdsname;
	auto lname = poVar->GetAttribute("long_name");
	if (lname) s.source_name_long = lname->ReadAsString();

	auto mval = poVar->GetAttribute("missing_value");
	if (mval) {
		s.m_missing_value = mval->ReadAsDouble();
		s.m_hasNA = true;
	} else {
		s.m_missing_value = poVar->GetNoDataValueAsDouble(&s.m_hasNA);
	}
	
	if (dimcount.size() < 2) {
		setError("insufficient dimensions");
		return false;
	}

	if (xyz.size() < 2) {
		xyz = {2,1,0};
	}

	int ix = ndim - 1;
	int iy = ix - 1;
	int it = 0;
	int iz = -1;
	if (ix == 3) {
		iz = 1;
		xyz = {ix, iy, iz, it};
	} else if (ndim == 2) {
		xyz = {ix, iy};
	}

	//Rcpp::Rcout << ix << ", " << iy << ", " << iz << ", " << it << std::endl;
	
	SpatExtent e;
 	s.ncol = dimcount[ix];
	s.nrow = dimcount[iy];
	s.nlyr = 1;

	// to do: check for equal spacing if x or y dim

	double ystart = dimvals[iy][0];
	double yend = dimvals[iy][dimvals[iy].size()-1];
	double res = (yend - ystart) / (s.nrow-1);
	e.ymax = yend + 0.5 * res;
	e.ymin = ystart - 0.5 * res;
	
	// reverse signaling to match non-md
	s.flipped = false;
	if (e.ymin > e.ymax) {
		std::swap(e.ymin, e.ymax);
		s.flipped = true;
	}
	s.m_dimnames.push_back(dimnames[iy]);

	double xstart = dimvals[ix][0];
	double xend = dimvals[ix][dimvals[ix].size()-1];
	res = (xend - xstart) / (s.ncol-1);
	e.xmin = dimvals[ix][0] - 0.5 * res;
	e.xmax = dimvals[ix][dimvals[ix].size()-1] + 0.5 * res;

	s.m_dimnames.push_back(dimnames[ix]);

	if (s.m_ndims > 2) {
		s.nlyr = dimcount[it];
		s.m_dimnames.push_back(dimnames[it]);
		std::string msg;
		parse_ncdf_time(s, dimunits[it], "standard", dimvals[it], msg);
	}

	if (iz >= 0) {
		s.m_dimnames.push_back(dimnames[iz]);
		s.depthname = dimnames[iz];
		s.depth = dimvals[iz];
		s.hasDepth = true;
		s.nlyr *= dimcount[iz];
	}
	
	for (size_t i=0; i<xyz.size(); i++) {
		s.m_dims.push_back(xyz[i]);
	}
	s.extent = e;

	s.nlyrfile = s.nlyr;
	s.resize(s.nlyr);
	s.layers.resize(s.nlyr);
    std::iota(s.layers.begin(), s.layers.end(), 0);
	s.rotated = false;
	s.memory = false;
	s.filename = fname;
	s.hasValues = true;
	s.unit = std::vector<std::string>(s.nlyr, poVar->GetUnit());
	s.is_multidim = true;

// layer names
	std::vector<std::string> nms;
	if (iz >= 0) {
		size_t niz = dimcount[iz];
		size_t ntm = dimcount[it];	
		nms.resize(ntm * niz, subdsname + "-");
		size_t k = 0;
		for (size_t i=0; i<niz; i++) {
			for (size_t j=0; j<ntm; j++) {
				nms[k] += std::to_string(j + 1) + "_" + s.depthname + "=" + std::to_string(dimvals[iz][i]);
				k++;
			}
		}
	} else {
		nms.resize(s.nlyr, subdsname + "-");
		for (size_t i=0; i<nms.size(); i++) {
			nms[i] += std::to_string(i + 1);
		}	
	}
	s.names = nms;

// time
	s.m_counts = dimcount;
	setSource(s);
	if (verbose) {
		for (size_t i=0; i<s.m_dims.size(); i++){
			Rcpp::Rcout << s.m_dims[i] << " " << dimnames[i] << " " << s.m_counts[i] << std::endl;
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

/*	
	Rcpp::Rcout << "dims: ";
	for (size_t i=0; i<dims.size(); i++) {Rcpp::Rcout << dims[i] << " ";}
	Rcpp::Rcout << "\n";

   std::iota(s.layers.begin(), s.layers.end(), 0);

*/

	offset[source[src].m_dims[0]] = col;
	offset[source[src].m_dims[1]] = row;
	size_t ndim = source[src].m_dims.size();
	std::vector<size_t> count(source[src].m_ndims, 1);
	count[source[src].m_dims[0]] = ncols;
	count[source[src].m_dims[1]] = nrows;

	std::vector<GPtrDiff_t> stride;
	if (!source[0].flipped) { 
		stride.resize(ndim, 1);
		stride[ndim-2] = -1;
		offset[ndim-2] = nrow() - row - 1;
	}
	
	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);
	if (source[src].in_order(false)) {
		if (ndim == 3) {
			offset[source[src].m_dims[2]] = source[0].layers[0];
			count[source[src].m_dims[2]] = source[0].layers.size();
		} else if (ndim == 4) {
			count[source[src].m_dims[2]] = source[0].depth.size();
			count[source[src].m_dims[3]] = source[0].time.size();		
		}
		size_t n=vprod(count, false);
		out.resize(n);
		GDALMDArrayRead(source[src].gdalmdarray, &offset[0], &count[0], &stride[0], NULL, hDT, &out[0], NULL, 0);
    } else {
		count[source[src].m_dims[2]] = 1;
		out.resize(0);
		out.reserve(ncols*nrows*source[0].layers.size());
		std::vector<double> lyr;
		size_t n=vprod(count, false);
		lyr.resize(n);
		for (size_t i=0; i<source[0].layers.size(); i++) {
			if (ndim == 3) {
				offset[source[src].m_dims[2]] = source[0].layers[i];
			} else if (ndim == 4) {
				setError("not handled yet");
				return false;
				count[source[src].m_dims[3]] = 1;		
			}
			GDALMDArrayRead(source[src].gdalmdarray, &offset[0], &count[0], NULL, NULL, hDT, &lyr[0], NULL, 0);
			out.insert(out.end(), lyr.begin(), lyr.end());
		}
	}
	GDALExtendedDataTypeRelease(hDT);


//	for (size_t i=0; i<offset.size(); i++) Rcpp::Rcout << offset[i] << ", "; 
//	Rcpp::Rcout << std::endl;
//	for (size_t i=0; i<count.size(); i++) Rcpp::Rcout << count[i] << ", "; 
//	Rcpp::Rcout << std::endl;


//	if (!source[0].flipped) { vflip(out, ncols); }

//	tpose(out, nrows, ncols, nlyr());

	if (source[src].m_hasNA) {
		Rcpp::Rcout << source[src].m_missing_value << std::endl;
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


/*
void tpose(std::vector<double> &v, size_t nr, size_t nc, size_t nl) {
	std::vector<double> vv(v.size());
	for (size_t lyr=0; lyr<nl; lyr++) {
		size_t off = lyr*nc*nr;
		for (size_t r = 0; r < nr; r++) {
			size_t rnc = off + r * nc;
			for (size_t c = 0; c < nc; c++) {
				vv[c*nr+r+off] = v[rnc+c];
			}
		}
	}
	v = vv;
}

void vflip(std::vector<double> &v, size_t nc) {
	std::vector<double> vv;
	vv.reserve(v.size());
	size_t nr = (v.size() / nc) - 1;
	for (int i=nr; i>=0; i--) {
		size_t b = i * nc;
		vv.insert(vv.end(), v.begin()+b, v.begin()+b+nc);
	}
	v = vv;
}

*/

