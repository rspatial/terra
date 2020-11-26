#include "spatRaster.h"
#include "time.h"
#include "recycle.h"

#include "string_utils.h"
#include "gdal_info.h"

bool good_ends(std::string const &s) {
	std::vector<std::string> end = {"_bnds", "_bounds", "lat", "lon", "longitude", "latitude"};
	for (size_t i=0; i<end.size(); i++) {
		if (s.length() >= end[i].length()) {
			if (s.compare(s.length() - end[i].length(), s.length(), end[i]) == 0) {
				return false;
			}
		} 
	}
	return true;
}

bool SpatRaster::constructFromNCDFsds(std::string filename, std::vector<std::string> meta, std::vector<int> subds, std::vector<std::string> subdsname) {

	std::vector<std::vector<std::string>> info = parse_metadata_sds(meta);
	int n = info[0].size();
	std::vector<std::string> sd, varname;

//	std::vector<unsigned> varnl;
// for selection based on nlyr
//	std::vector<size_t> nl(n);
//	for (size_t i=0; i<nl.size(); i++) {
//		nl[i] = stol(info[5][i]);
//	}
	
	if (info[0].size() == 0) {
		return false;
	}
	// select sds by index
	if (subds[0] >=0) {
		for (size_t i=0; i<subds.size(); i++) {
			if (subds[i] >=0 && subds[i] < n) {
				sd.push_back(info[0][subds[i]]);
				varname.push_back(info[1][i]);
				//varnl.push_back(nl[i]);			
			} else {
				std::string emsg = std::to_string(subds[i]+1) + " is not valid. There are " + std::to_string(info[0].size()) + " subdatasets\n";
				setError(emsg);
				return false;
			}
		}
	// select by name	
	} else if (subdsname[0] != "") {
		for (size_t i=0; i<subdsname.size(); i++) {
			int w = where_in_vector(subdsname[i], info[1]);
			if (w >= 0) {
				sd.push_back(info[0][w]);
				varname.push_back(info[1][w]);
				//varnl.push_back(nl[w]);			
			} else {
				std::string emsg = concatenate(info[1], ", ");
				emsg = subdsname[i] + " not found. Choose one of:\n" + emsg;
				setError(emsg);
				return false;
			}
		}
	// select all
	} else {
		// eliminate sources based on names like "*_bnds" and "lat"
		for (size_t i=0; i<info[1].size(); i++) {
			if (good_ends(info[1][i])) {
				sd.push_back(info[0][i]);
				varname.push_back(info[1][i]);
			} 
		}
	// perhaps also remove sds with 1 row or col of other sds have more.		
/*
// alternatively, eliminate based on nlyr()
		size_t mxnl = *max_element(nl.begin(), nl.end());
		for (size_t i=0; i<nl.size(); i++) {
			if (nl[i] == mxnl) {
				sd.push_back(info[0][i]);
				varname.push_back(info[1][i]);
				varnl.push_back(nl[i]);			
			}			
		}
*/		
	}
	
	std::vector<size_t> srcnl;
	bool success = constructFromFile(sd[0], {-1}, {""});
	if (!success) {
		return false;
	}
	std::vector<std::string> skipped, used;
	srcnl.push_back(nlyr());
	used.push_back(varname[0]);				
	SpatRaster out;
    for (size_t i=1; i < sd.size(); i++) {
//		printf( "%s\n", sd[i].c_str() );
		success = out.constructFromFile(sd[i], {-1}, {""});
		if (success) {
			if (out.compare_geom(*this, false, false)) {
				addSource(out);
				srcnl.push_back(out.nlyr());
				used.push_back(varname[i]);				
			} else {
				skipped.push_back(varname[i]);
			}
		} else {
			skipped.push_back(varname[i]);
		}
	}

	if (skipped.size() > 0) {
		std::string s="skipped sub-datasets (different geometry, see 'describe_sds'):\n";
		for (size_t i=0; i<skipped.size(); i++) {
			s += skipped[i];
			if ((i%3) == 0) {
				s += "\n";
			} else {
				if (i > 0) s += ", ";
			}
		}
		addWarning(s);
	}


	std::vector<std::string> lyrnames;
	for (size_t i=0; i<used.size(); i++) {
		std::vector<std::string> nms = {used[i]};
		recycle(nms, srcnl[i]);
		make_unique_names(nms);
		lyrnames.insert(lyrnames.end(), nms.begin(), nms.end());
		//Rcpp::Rcout << used[i] << std::endl;
		//Rcpp::Rcout << nms.size() << std::endl;
		
	}
	if (lyrnames.size() > 0) {
		success = setNames(lyrnames, false);
	}

	return true;
}



std::vector<int_64> str2int64v(std::string s, std::string delim) {
	std::vector<int_64> out;
	size_t pos = 0;
	while ((pos = s.find(delim)) != std::string::npos) {
		std::string v = s.substr(0, pos);
		s.erase(0, pos + 1);
		out.push_back(std::stoll(v));
	}
	out.push_back(std::stoll(s));
	return out;
}


std::vector<int_64> ncdf_time(const std::vector<std::string> &metadata) {

	bool fu=false;
	bool fv=false;
	bool fc=false;
	std::string origin, values, calendar;

	for (size_t i=0; i<metadata.size(); i++) {
		if (!fc) {
			std::string pattern = "time#calendar=";
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				calendar = metadata[i];
				calendar.erase(calendar.begin(), calendar.begin()+pattern.size());  
				fc = true;
			}
		}
		if (!fu) {
			std::string pattern = "time#units=";
			std::size_t found = metadata[i].find(pattern);
			if (found != std::string::npos) {
				origin = metadata[i];
				origin.erase(origin.begin(), origin.begin()+pattern.size());  
				fu = true;
			}
		}
		if (!fv) {
			std::size_t found = metadata[i].find("NETCDF_DIM_time_VALUES=");
			if (found != std::string::npos) {
				values = metadata[i].substr(24, metadata[i].size()-1);  
				if (values.size() == 0) {
					break;
				}
				fv = true;
			}
		}
		if (fc & fu & fv) break;
	}

	std::vector<int_64> out, bad;
	if (!(fu & fv)) {
		return out;
	}


	bool days=false; 
	bool hours=false;
	//bool seconds = false; 

	out = str2int64v(values, ",");

	if ((origin.find("seconds")) != std::string::npos) {
		//seconds = true;
	} else if ((origin.find("hours")) != std::string::npos) {
		hours = true;
	} else if ((origin.find("days")) != std::string::npos) {
		days = true;
	}

			
	//if (calendar =="gregorian" || calendar =="proleptic_gregorian" || calendar=="standard" || calendar == "") { // ok }
	//if ((calendar == "noleap") || (calendar == "365 day") || (calendar == "365_day") || (calendar == "360 day") || (calendar == "360_day")) {
	
	bool found = false;
	size_t pos;
	if ((pos = origin.find("from")) != std::string::npos) {
		origin.erase(0, pos + 5);
		found = true;
	} else if ((pos = origin.find("since")) != std::string::npos) {
		origin.erase(0, pos + 6);
		found = true;	
	}
	SpatTime_t offset = 0;
	if (found) {
		if (days) {
			if (calendar == "noleap" || calendar == "365_day" || calendar == "365 day") { 
				std::vector<int> ymd = getymd(origin);
				for (int_64 &d : out) d = time_from_day_noleap(ymd[0], ymd[1], ymd[2], d);
			} else if (calendar == "360_day" || calendar == "360 day") { 
				std::vector<int> ymd = getymd(origin);
				for (int_64 &d : out) d = time_from_day_360(ymd[0], ymd[1], ymd[2], d);
			} else {
				std::vector<int> ymd = getymd(origin);
				for (int_64 &d : out) d = time_from_day(ymd[0], ymd[1], ymd[2], d);
			}
		} else if (hours) {
			hours_to_time(out, origin);

			//std::vector<int> ymd = getymd(origin);
			//for (int_64 &d : out) d = time_from_hour(ymd[0], ymd[1], ymd[2], d);
		} else { // seconds
			offset = get_time_string(origin);
			for (int_64 &d : out) d = d + offset;
		}		
	}

	return out;
}
	

void RasterSource::set_names_time_ncdf(std::vector<std::string> metadata, std::vector<std::vector<std::string>> bandmeta) {

	//std::vector<std::string> varname getmeta("NETCDF_VARNAME");
	//NETCDF_DIM_time=729649.5
	//NETCDF_VARNAME=pr
	
	//varname = ncdf_name(metadata);
		
	std::vector<int_64> x = ncdf_time(metadata);
	if (x.size() == nlyr) {
		time = x;
		hasTime = true;
	}
}

	

/*

std::string ncdf_name(const std::vector<std::string> &metadata) {
	bool wasfound = false;
	std::string name = "";
	for (size_t i=0; i<metadata.size(); i++) {
		std::size_t found = metadata[i].find("NETCDF_DIM_");
		if (found == std::string::npos) {
			if (wasfound) {
				std::size_t pos = metadata[i].find("#");
				name = metadata[i].substr(0, pos);
				break;
			}
		} else {
			wasfound = true;
		}
	}
	return name;
}


std::vector<std::vector<std::string>> metatime(std::vector<std::string> meta) {
	std::vector<std::vector<std::string>> out(meta.size());
	std::string delim = "=";
	for (size_t i=0; i<meta.size(); i++) {
		std::string s = meta[i];
		size_t pos = s.find(delim);
		if (pos != std::string::npos) {
			out[i].push_back(s.erase(pos+1, std::string::npos));
			out[i].push_back(s.erase(0, pos+1));
		} else {
			out[i].push_back(s);
		}
	}
	return out;
}


*/

