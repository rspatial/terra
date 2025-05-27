
#include "spatRaster.h"


#if GDAL_VERSION_MAJOR >= 3 && GDAL_VERSION_MINOR >= 4

#include "proj.h"

#include "ogr_spatialref.h"

#include "gdal_priv.h"
#include "gdal.h"
#include "crs.h"
#include "string_utils.h"
#include "vecmath.h"


std::vector<std::string> GetArrayNames(std::shared_ptr<GDALGroup> x) {
// FROM GDAL 3.11 (while not widely available).
// * Author:   Even Rouault <even.rouault at spatialys.com>
// * Copyright (c) 2019, Even Rouault <even.rouault at spatialys.com>
 
    std::vector<std::string> ret;
    std::list<std::shared_ptr<GDALGroup>> stackGroups;
    stackGroups.push_back(nullptr);  // nullptr means this
    while (!stackGroups.empty()) {
        std::shared_ptr<GDALGroup> groupPtr = std::move(stackGroups.front());
        stackGroups.erase(stackGroups.begin());
        const GDALGroup *poCurGroup = groupPtr ? groupPtr.get() : x.get();
        for (const std::string &arrayName :
             poCurGroup->GetMDArrayNames(nullptr))
        {
            std::string osFullName = poCurGroup->GetFullName();
            if (!osFullName.empty() && osFullName.back() != '/')
                osFullName += '/';
            osFullName += arrayName;
            ret.push_back(std::move(osFullName));
        }
        auto insertionPoint = stackGroups.begin();
        for (const auto &osSubGroup :
             poCurGroup->GetGroupNames(nullptr))
        {
            auto poSubGroup = poCurGroup->OpenGroup(osSubGroup);
            if (poSubGroup)
                stackGroups.insert(insertionPoint, std::move(poSubGroup));
        }
    }

    return ret;
}


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


std::vector<std::string> ncdf_keep(std::vector<std::string> const &s) {
	std::vector<std::string> out;
	out.reserve(s.size());
	std::vector<std::string> end = {"_bnds", "_bounds", "lat", "lon", "longitude", "latitude"};
	for (size_t j=0; j<s.size(); j++) {
		bool add = true;
		for (size_t i=0; i<end.size(); i++) {
			if (s[j].length() >= end[i].length()) {
				if (s[j].compare(s[j].length() - end[i].length(), s[j].length(), end[i]) == 0) {
					add = false;
					continue;
				}
			}
		}
		if (add && (!(s[j] == "/x" || s[j] == "/y" || s[j] == "/northing" || s[j] == "/easting" || s[j] == "/time"))) {
			out.push_back(s[j]);
		}
	}
	return out;
}


void prints(std::vector<std::string> &x) {
	for (size_t i=0; i<x.size(); i++) {Rcpp::Rcout << x[i] << " ";}
	Rcpp::Rcout << "\n";	
}


bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> dims, bool noflip, bool guessCRS, std::vector<std::string> domains) {

	SpatRasterSource s;

	bool verbose = false;

    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(fname.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !poDataset ) {
		setError("not a good dataset");
        return false;
    }

/*
	if (subname.size() > 0) {
		std::string s = subname[0];
		if (s[0] == "/") s = s.substr(1,  s.length());
		if (in_string(s, "/") {
			std::vector<std::string> ss = strsplit_last(s, "/");
			subname[0] = ss[1];
			group = ss[0];
		}
	}
*/

	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no root group");
		return false;
    }
	
		
//	std::vector<std::string> x = poRootGroup->GetGroupNames();
//	prints(x);

// GDAL 3.11		
//	x = poRootGroup->GetMDArrayFullNamesRecursive();
//	prints(x);

/*
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
//				if (verbose) {
					Rcpp::Rcout << names[i] << ": ";	
					for (size_t j=0; j<ni; j++) {
						Rcpp::Rcout << dim_names[i][j] << " (" << dim_size[i][j] << ") ";	
					}
					Rcpp::Rcout << std::endl;
//				}
			}
		}
	}
*/

//	if (xyz.size() != 3) {
//		setError("you must supply three dimension indices");
//       return false;
//	}


	s.m_arrayname = "";
	if (!subname[0].empty()) {
		s.m_arrayname = subname[0];
		//int w = where_in_vector(s.m_arrayname, ar_names, false);
		//if (w < 0) {
			//setError("array " + s.m_arrayname + " not found. Should be one of:\n  " + concatenate(ar_names, ", "));
			//return false;
		//} 
//	} else if (sub.size() > 0) {
		//if ((sub[0] < 0) || (sub[0] >= (int)ar_names.size())) {
		//	setError("array number is out or range");
		//	return false;
		//} else {
		//	s.m_arrayname = ar_names[sub[0]];
		//}
	} else {
		std::vector<std::string> anms = ncdf_keep(GetArrayNames(poRootGroup));
		//prints(anms);
		s.m_arrayname = anms[anms.size()-1];
		if (anms.size() > 1)  {
			anms = {anms.begin(), anms.end() - 1}; 
			addWarning("using: " + s.m_arrayname + ". Other arrays are: \n" + concatenate(anms, "\n"));
		}
	}

	std::string startgroup="";
	auto poVar = poRootGroup->ResolveMDArray(s.m_arrayname.c_str(), startgroup, nullptr);
//    auto poVar = poRootGroup->OpenMDArray(s.m_arrayname.c_str());
    if( !poVar )   {
		setError("cannot find: " + s.m_arrayname);
		return false;
    }
	s.m_arrayname = poVar->GetFullName();

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
	if (guessCRS && wkt.empty()) {
		if (s.extent.xmin >= -181 && s.extent.xmax <= 361 && s.extent.ymin >= -91 && s.extent.ymax <= 91) {
			wkt = "OGC:CRS84";
			s.parameters_changed = true;
		}
	}
	std::string msg = "";
	if (!s.srs.set({wkt}, msg)) {
		addWarning(msg);
	}

// dimensions 
	std::vector<size_t> dimcount;
	std::vector<std::string> dimnames, dimunits;
	std::vector<std::vector<double>> dimvals;
	dimvals.reserve(4);

	std::string calendar = "";	
	std::vector<std::shared_ptr<GDALDimension>> dimData = poVar->GetDimensions();
	size_t ndim = dimData.size();
    for (size_t i=0; i<ndim; i++) {
		size_t n = dimData[i]->GetSize();
        dimcount.push_back(n);
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
	if (s.m_ndims < 2) {
		setError("insufficient dimensions");
		return false;
	}

	s.source_name = s.m_arrayname;
	auto lname = poVar->GetAttribute("long_name");
	if (lname) s.source_name_long = lname->ReadAsString();

	auto mval = poVar->GetAttribute("missing_value");
	if (mval) {
		s.m_missing_value = mval->ReadAsDouble();
		s.m_hasNA = true;
	} else {
		s.m_missing_value = poVar->GetNoDataValueAsDouble(&s.m_hasNA);
	}
	

	if (dims.size() < 2) {
		dims = {2,1,0};
	}

	int ix = ndim - 1;
	int iy = ix - 1;
	int it = 0;
	int iz = -1;
	if (ix == 3) {
		iz = 1;
		dims = {ix, iy, iz, it};
	} else if (ndim == 2) {
		dims = {ix, iy};
	}

	//Rcpp::Rcout << ix << ", " << iy << ", " << iz << ", " << it << std::endl;
	
	SpatExtent e;
 	s.ncol = dimcount[ix];
	s.nrow = dimcount[iy];
	s.nlyr = 1;

	// to do: check for equal spacing if x or y dim

	double start = dimvals[ix][0];
	double end = dimvals[ix][dimvals[ix].size()-1];
	double res = (end - start) / (s.ncol-1);
	e.xmin = start - 0.5 * res;
	e.xmax = end + 0.5 * res;

	start = dimvals[iy][0];
	end = dimvals[iy][dimvals[iy].size()-1];
	res = (end - start) / (s.nrow-1);
	e.ymax = end + 0.5 * res;
	e.ymin = start - 0.5 * res;
	
	s.flipped = false;
	if (e.ymin > e.ymax) {
		std::swap(e.ymin, e.ymax);
		s.flipped = true;
	}
//	s.m_names.push_back(dimnames[ix]);
//	s.m_names.push_back(dimnames[iy]);


	if (s.m_ndims > 2) {
		s.nlyr = dimcount[it];
//		s.m_names.push_back(dimnames[it]);
		std::string msg;
		parse_ncdf_time(s, dimunits[it], "standard", dimvals[it], msg);
	}

	if (iz >= 0) {
//		s.m_names.push_back(dimnames[iz]);
		s.depthname = dimnames[iz];
		s.depth = dimvals[iz];
		s.hasDepth = true;
		s.nlyr *= dimcount[iz];
	}
	s.nlyrfile = s.nlyr;
	s.resize(s.nlyr);
	s.layers.resize(s.nlyr);
    std::iota(s.layers.begin(), s.layers.end(), 0);
	
	for (size_t i=0; i<dims.size(); i++) {
		s.m_dims.push_back(dims[i]);
	}
	s.extent = e;

	bool app_so = true;
	size_t opsz = options.size();
	if (opsz > 0) {
		if (options[opsz-1] == "so=false") {
			app_so = false;
			options.resize(opsz-1); 
		}
	}

	if (app_so) {
		bool hasScale=false;
		double scale = poVar->GetScale(&hasScale, nullptr);
		if (scale == 1) hasScale = false;
		bool hasOffset=false;
		double offset = poVar->GetOffset(&hasOffset, nullptr);
		if (offset == 0) hasOffset = false;
		if (hasScale || hasOffset) {
			s.has_scale_offset = std::vector<bool>(s.nlyr, true);
			s.offset = std::vector<double>(nlyr(), offset);
			s.scale = std::vector<double>(nlyr(), scale);
		}
	}
	
	s.rotated = false;
	s.memory = false;
	s.filename = fname;
	s.hasValues = true;
	s.unit = std::vector<std::string>(s.nlyr, poVar->GetUnit());
	s.is_multidim = true;

// layer names
	std::vector<std::string> nms;
	std::vector<std::string> arn = strsplit_last(s.m_arrayname, "/");
	std::string arname = arn[arn.size()-1];
	if (iz >= 0) {
		size_t niz = dimcount[iz];
		size_t ntm = dimcount[it];	
		nms.resize(ntm * niz, arname + "-");
		size_t k = 0;
		for (size_t i=0; i<niz; i++) {
			std::string sz = double_to_string(dimvals[iz][i]);
			for (size_t j=0; j<ntm; j++) {
				nms[k] += std::to_string(j + 1) + "_" + s.depthname + "=" + sz;
				k++;
			}
		}
	} else {
		nms.resize(s.nlyr, arname + "-");
		for (size_t i=0; i<nms.size(); i++) {
			nms[i] += std::to_string(i + 1);
		}	
	}
	s.names = nms;

// time
	s.m_size = dimcount;
	s.m_names = dimnames;
	
	setSource(s);
	if (verbose) {
		for (size_t i=0; i<s.m_dims.size(); i++){
			Rcpp::Rcout << s.m_dims[i] << " " << dimnames[i] << " " << s.m_size[i] << std::endl;
		}
	}
	return true;
}



bool SpatRaster::readStartMulti(size_t src) {

//	Rcpp::Rcout << "readStartMulti\n";
/*
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

	GDALMDArrayH hVar = GDALGroupOpenMDArray(hGroup, source[src].m_arrayname.c_str(), NULL);
    GDALGroupRelease(hGroup);
    if (!hVar) {
		setError("array '" + source[src].m_arrayname + "' is not available");
		return false;
    }
*/

    auto poDataset = std::unique_ptr<GDALDataset>(GDALDataset::Open(source[src].filename.c_str(), GDAL_OF_MULTIDIM_RASTER ));
    if( !poDataset ) {
		setError("not a good dataset");
        return false;
    }


	std::shared_ptr<GDALGroup> poRootGroup = poDataset->GetRootGroup();
    if( !poRootGroup ) {
		setError("no roots");
		return false;
    }
//    GDALReleaseDataset(hDS);

//	std::string startgroup="";
//	auto poVar = poRootGroup->ResolveMDArray(source[src].m_arrayname.c_str(), startgroup, nullptr);

    auto poVar = poRootGroup->OpenMDArray(source[src].m_arrayname.c_str());
    if( !poVar )   {
		setError("cannot find: " + source[src].m_arrayname);
		return false;
    }


	if (source[src].has_scale_offset[0]) {
		source[src].m_array = poVar->GetUnscaled();
	} else {
		source[src].m_array = poVar;
	}
	source[src].open_read = true;
	return true;
}


bool SpatRaster::readStopMulti(size_t src) {
//	Rcpp::Rcout << "readStopMulti\n";
	source[src].open_read = false;
	return true;
}



bool SpatRaster::readChunkMulti(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {

//	Rcpp::Rcout << "readChunkMulti\n";
	std::vector<GUInt64> offset(source[src].m_ndims, 0);

/*	
	std::vector<size_t> dims = source[src].m_dims;
	Rcpp::Rcout << "dims: ";
	for (size_t i=0; i<dims.size(); i++) {Rcpp::Rcout << dims[i] << " ";}
	Rcpp::Rcout << "\n";
*/

	offset[source[src].m_dims[0]] = col;
	offset[source[src].m_dims[1]] = row;
	size_t ndim = source[src].m_dims.size();
	std::vector<size_t> count(source[src].m_ndims, 1);
	count[source[src].m_dims[0]] = ncols;
	count[source[src].m_dims[1]] = nrows;

	std::vector<GPtrDiff_t> stride;
	if (!source[src].flipped) { 
		stride.resize(ndim, 1);
		stride[ndim-2] = -1;
		offset[ndim-2] = nrow() - row - 1;
	}

//	GDALExtendedDataTypeH hDT = GDALExtendedDataTypeCreate(GDT_Float64);
//	std::vector<double> out;
	size_t insize = data.size();

	auto dt = GDALExtendedDataType::Create(GDT_Float64);

	if (source[src].in_order(false)) {
		if (ndim == 3) {
			offset[source[src].m_dims[2]] = source[src].layers[0];
			count[source[src].m_dims[2]] = source[src].layers.size();
		} else if (ndim == 4) {
			count[source[src].m_dims[2]] = source[src].depth.size();
			count[source[src].m_dims[3]] = source[src].time.size();		
		}
		size_t n=vprod(count, false);
		data.resize(insize + n);
		source[src].m_array->Read(&offset[0], &count[0], &stride[0], NULL, dt, &data[insize], NULL, 0);
    } else {
		count[source[src].m_dims[2]] = 1;
//		out.resize(0);
//		out.reserve(ncols*nrows*source[src].layers.size());

		data.resize(insize + ncols*nrows*source[src].layers.size());
//		size_t n=vprod(count, false);
		for (size_t i=0; i<source[src].layers.size(); i++) {
			if (ndim == 3) {
				offset[source[src].m_dims[2]] = source[src].layers[i];
			} else if (ndim == 4) {
				setError("not handled yet");
				return false;
				count[source[src].m_dims[3]] = 1;		
			}
			source[src].m_array->Read(&offset[0], &count[0], NULL, NULL, dt, &data[insize], NULL, 0);
		}
	}
//	GDALExtendedDataTypeRelease(hDT);

//	for (size_t i=0; i<offset.size(); i++) Rcpp::Rcout << offset[i] << ", "; 
//	Rcpp::Rcout << std::endl;
//	for (size_t i=0; i<count.size(); i++) Rcpp::Rcout << count[i] << ", "; 
//	Rcpp::Rcout << std::endl;


//	if (!source[0].flipped) { vflip(out, ncols); }

//	tpose(out, nrows, ncols, nlyr());

	if (source[src].m_hasNA) {
//		Rcpp::Rcout << source[src].m_missing_value << std::endl;
		std::replace (data.begin()+insize, data.end(), source[src].m_missing_value, (double)NAN);
	}

//	data.insert(data.end(), out.begin(), out.end());
	return true;
}

//#include "gdalio.h"

bool SpatRaster::writeStartMulti(SpatOptions &opt, const std::vector<std::string> &srcnames) {

	if (!hasValues()) {
		setError("there are no cell values");
		return false;
	}

	std::string filename = opt.get_filename();	
	if (filename.empty()) {
		setError("empty filename");
		return(false);
	} 
	// assure filename won't be used again
	opt.set_filenames({""});
	std::string driver = "netCDF";
/*
	std::string driver = opt.get_filetype();
	getGDALdriver(filename, driver);
	if (driver.empty()) {
		setError("cannot guess file type from filename");
		return(false);
	}
	if (driver != "netCDF") {
		if (driver.empty()) {
			setError("multi-dim only implemented for netCDF");
			return(false);
		}
	}
*/
	
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName(driver.c_str());

	GDALDataset *poDS;
    poDS = poDriver->CreateMultiDimensional(filename.c_str(), NULL, NULL);
	if (poDS == NULL) {
		setError("failed writing "+ driver + " file");
		GDALClose( (GDALDatasetH) poDS );
		return false;
	}

    auto rg = poDS->GetRootGroup();
	
	std::vector<std::shared_ptr<GDALDimension>> dim_ptrs;
	auto dt = GDALExtendedDataType::Create(GDT_Float64);


	size_t nz = nlyr();
	size_t ny = nrow();
	size_t nx = ncol();
	dim_ptrs.push_back(rg->CreateDimension("Z", "", "", nz));
	dim_ptrs.push_back(rg->CreateDimension("Y", "", "", ny));
	dim_ptrs.push_back(rg->CreateDimension("X", "", "", nx));


    auto var = rg->CreateMDArray("Z", {dim_ptrs[0]}, dt);
    var = rg->OpenMDArray("Z");

	std::vector<double> dvals(nz);
	std::iota(dvals.begin(), dvals.end(), 0);

	std::vector<size_t> count = {nz};
	std::vector<GUInt64> start = {0};

	var->Write(start.data(), count.data(), nullptr, nullptr, dt, &dvals[0]); 

    var = rg->CreateMDArray("Y", {dim_ptrs[1]}, dt);
	yFromRow(dvals);
	count = {ny};
	var->Write(start.data(), count.data(), nullptr, nullptr, dt, &dvals[0]); 

    var = rg->CreateMDArray("X", {dim_ptrs[2]}, dt);
	xFromCol(dvals);
	count = {nx};
	var->Write(start.data(), count.data(), nullptr, nullptr, dt, &dvals[0]); 
	
	std::string vname = source[0].source_name.empty() ? "array" : source[0].source_name;
	
    var = rg->CreateMDArray(vname, dim_ptrs, GDALExtendedDataType::Create(GDT_Float64));
	

	std::string wkt = source[0].srs.wkt;
	if (!wkt.empty()) {
		OGRSpatialReference *srs = NULL;
		srs = new OGRSpatialReference;
		const char *cp = wkt.c_str();
		srs->importFromWkt(cp);
		if (srs != NULL) {
			if (!var->SetSpatialRef(srs)) {
				addWarning("failed to assign CRS to array");
			}
			CPLFree(srs);
		}
	}

	source[0].m_array = var;
	source[0].gdalconnection = poDS;

	return true;
}

bool SpatRaster::writeValuesMulti(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols){
	SpatOptions opt;
	std::vector<GUInt64> start = {0, 0, 0};
	std::vector<size_t> count = {nlyr(), nrow(), ncol()};
	source[0].m_array->Write(start.data(), count.data(), nullptr, nullptr, GDALExtendedDataType::Create(GDT_Float64), &vals[0]); 
	return true;
}

bool SpatRaster::writeStopMulti() {
	//GDALMDArrayRelease(
	source[0].m_array.reset();
	GDALClose( source[0].gdalconnection );
	return true;
}

#else


bool SpatRaster::constructFromFileMulti(std::string fname, std::vector<std::string> subname, std::vector<std::string> drivers, std::vector<std::string> options, std::vector<int> dims, bool noflip, bool guessCRS, std::vector<std::string> domains) {
	setError("multidim is not supported with GDAL < 3.4");
	return false;
}

bool SpatRaster::readStartMulti(size_t src) {
	return false;
}

bool SpatRaster::readStopMulti(size_t src) {
	return false;
}

bool SpatRaster::readChunkMulti(std::vector<double> &data, size_t src, size_t row, size_t nrows, size_t col, size_t ncols) {
	return false
}

bool SpatRaster::writeStartMulti(SpatOptions &opt, const std::vector<std::string> &srcnames) {
	return false
}

bool SpatRaster::writeValuesMulti(std::vector<double> &vals, size_t startrow, size_t nrows, size_t startcol, size_t ncols){
	return false;
}

bool SpatRaster::writeStopMulti() {
	return false
}

#endif


std::vector<double> SpatRaster::readValuesMulti(size_t src, size_t row, size_t nrows, size_t col, size_t ncols, int lyr) {

	Rcpp::Rcout << "readValuesMulti\n";
	std::vector<double> out;
	if (lyr < 0) {
		if (!readStartMulti(src)) {
			return out;
		}
		readChunkMulti(out, src, row, nrows, col, ncols);
		readStopMulti(src);
		return out;
	} else {
		Rcpp::Rcout << "empty\n";
		return out;
	}
}


bool SpatRaster::readRowColMulti(size_t src, std::vector<std::vector<double>> &out, size_t outstart, std::vector<int_64> &rows, const std::vector<int_64> &cols) {

//	Rcpp::Rcout << "readRowColMulti " << src << "\n";
	if (!readStartMulti(src)) {
		return false;
	}
	size_t n = rows.size();
	size_t nl = source[src].layers.size();
	size_t outend = outstart + nl;
	for (size_t i=outstart; i<outend; i++) {
		out[i].reserve(n); // = std::vector<double> (n, NAN);
	}
	
	std::vector<double> value;
	for (size_t i=0; i<n; i++) {
//		Rcpp::Rcout << rows[i] << " " << cols[i] << " ";
		if ((rows[i] < 0) || (cols[i] < 0)) {
			value.resize(value.size() + nl, NAN);
		} else {
			if (!readChunkMulti(value, src, rows[i], 1, cols[i], 1)) {
				return false;
			}
		}
	}
	readStopMulti(src);
	for (size_t i=0; i<n; i++) {
		out[outstart+(i%nl)].push_back(value[i]);
	}
	return true;
}


void getSampleRowCol2(std::vector<int_64> &oldrow, std::vector<int_64> &oldcol, size_t nrows, size_t ncols, size_t snrow, size_t sncol) {

	double rf = nrows / (double)(snrow);
	double cf = ncols / (double)(sncol);
	//double rstart = std::floor(0.5 * rf);
	//double cstart = std::floor(0.5 * cf);
	double rstart = 0.5 * rf;
	double cstart = 0.5 * cf;
	
	std::vector<int_64> xcol, xrow;
	xcol.reserve(sncol);
	for (size_t i =0; i<sncol; i++) {
        xcol.push_back(i * cf + cstart);
	}
	xrow.reserve(snrow);
	for (size_t i =0; i<snrow; i++) {
        xrow.push_back(i * rf + rstart);
	}
	oldrow.reserve(sncol * snrow);
	oldcol.reserve(sncol * snrow);
	for (size_t i =0; i<snrow; i++) {
		for (size_t j=0; j<sncol; j++) {
			oldrow.push_back(xrow[i]);
			oldcol.push_back(xcol[j]);
		}
	}
}


std::vector<double> SpatRaster::readSampleMulti(size_t src, size_t srows, size_t scols, bool overview) {
//	Rcpp::Rcout << "readSampleMulti\n";
	std::vector<int_64> colnr, rownr;
	getSampleRowCol2(rownr, colnr, nrow(), ncol(), srows, scols);
	std::vector<std::vector<double>> out(source[src].layers.size());
	readRowColMulti(src, out, 0, rownr, colnr);
    return out[0];
}


SpatRaster SpatRaster::writeRasterM(SpatOptions &opt) {
	SpatRaster out;
	
	std::vector<std::string> fnames = opt.get_filenames();

	if (!writeStartMulti(opt, {""})) {
		out.setError(getError());
		return out;
	}
	std::vector<double> vals = getValues(-1, opt);
	if (!writeValuesMulti(vals, 0, nrow(), 0, ncol())) {
		out.setError(getError());
		return out;		
	}
	writeStopMulti();

	std::vector<std::string> empty;
	std::vector<int> dims = {-1};
	
	out.constructFromFileMulti(fnames[0], empty, empty, empty, dims, false, false, {""});
	return out;
}


std::vector<bool> SpatRaster::is_multidim() {
	std::vector<bool> out;
	out.reserve(source.size());
	for (size_t i=0; i<source.size(); i++) {
		out.push_back(source[i].is_multidim);
	}
	return(out);
}

std::vector<std::vector<std::string>> SpatRaster::dim_names() {
	std::vector<std::vector<std::string>> out(source.size());
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].is_multidim) {
			out[i] = source[i].m_names;
		}
	}
	return(out);
	
}

std::vector<std::vector<size_t>> SpatRaster::dim_order() {
	std::vector<std::vector<size_t>> out(source.size());
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].is_multidim) {
			out[i] = source[i].m_order;
		}
	}
	return out;
}

std::vector<std::vector<size_t>> SpatRaster::dim_size() {
	std::vector<std::vector<size_t>> out(source.size());
	for (size_t i=0; i<source.size(); i++) {
		if (source[i].is_multidim) {
			out[i] = source[i].m_size;
		}
	}
	return out;
}

