// Copyright (c) 2018-2020  Robert J. Hijmans
//
// This file is part of the "spat" library.
//
// spat is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// spat is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with spat. If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
//#include "spatMessages.h"

#ifndef standalone
	#define useRcpp
#endif


#ifndef nogdal
  #define useGDAL
#endif


/*
#ifdef useGDAL
	#include "gdal_priv.h"
#endif
*/

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif




class SpatMessages {
	public:
		bool has_error = false;
		bool has_warning = false;
		std::string error;
		std::vector<std::string> warnings;

		void setError(std::string s) {
			has_error = true;
			error = s;
		}

		std::string getError() {
			has_error = false;
			return error;
		}
		
		void addWarning(std::string s) {
			has_warning = true;
			warnings.push_back(s);
		}

		std::string getWarnings() {
			std::string w = "";
			for (size_t i = 0; i<warnings.size(); i++) {
				w += warnings[i] + "\n" ;
			}
			warnings.resize(0);
			has_warning = false;
			return w;
		}
		
		std::vector<std::string> getMessages() {
			std::string warns = getWarnings();
			std::string error = getError();
			std::vector<std::string> msg = { error, warns};
			return msg;
		}
};


class SpatOptions {
	private:
		std::string tempdir = "";
		bool todisk = false;
		double memfrac = 0.6;

	public:
		std::string def_datatype = "FLT4S";
		std::string def_filetype = "GTiff";
		//std::string def_bandorder = "BIL";
		bool overwrite = false;
		unsigned progress = 3;
		unsigned blocksizemp = 4;
		size_t steps = 0;
		double NAflag = NAN;
		bool def_verbose = false;
		bool verbose = false;
		std::string datatype = "";
		//std::string bandorder = "";
		std::string filetype = "";
		std::vector<std::string> filenames = {""};
		std::vector<std::string> gdal_options;
		std::vector<std::string> names;

		SpatOptions();
		SpatOptions(const SpatOptions &opt);
		SpatOptions deepCopy(const SpatOptions &opt);

		// permanent
		bool get_todisk();
		void set_todisk(bool b);
		double get_memfrac();
		void set_memfrac(double d);
		std::string get_tempdir();
		void set_tempdir(std::string d);

		std::string get_def_datatype();
		std::string get_def_bandorder();
		std::string get_def_filetype();
		bool get_def_verbose();
		void set_def_datatype(std::string d);
		//void set_def_bandorder(std::string d);
		void set_def_filetype(std::string d);

		// single use
		void set_verbose(bool v);
		void set_def_verbose(bool v);
		void set_NAflag(double flag);
		//void set_filename(std::string f);
		void set_filenames(std::vector<std::string> f);
		void set_filetype(std::string d);
		void set_datatype(std::string d);
		//void set_bandorder(std::string d);
		void set_overwrite(bool b);
		void set_progress(unsigned p);
		void set_blocksizemp(unsigned x);
		std::string get_filename();
		std::vector<std::string> get_filenames();
		std::string get_filetype();
		std::string get_datatype();
		//std::string get_bandorder();
		bool get_verbose();
		double get_NAflag();
		bool get_overwrite();
		unsigned get_progress();
		bool do_progress(unsigned n);
		unsigned get_blocksizemp();
		void set_steps(size_t n);
		size_t get_steps();

		SpatMessages msg;
};



class SpatExtent {
	public:
		double xmin, xmax, ymin, ymax;
		double inf = std::numeric_limits<double>::infinity();
		double neginf = -std::numeric_limits<double>::infinity();
//		SpatExtent() {xmin = inf; xmax = neginf; ymin = inf; ymax = neginf;}
		SpatExtent() {xmin = -180; xmax = 180; ymin = -90; ymax = 90;}
		SpatExtent(double _xmin, double _xmax, double _ymin, double _ymax) {xmin = _xmin; xmax = _xmax; ymin = _ymin; ymax = _ymax;}

		void intersect(SpatExtent e) { // check first if intersects?
			xmin = std::max(xmin, e.xmin);
			xmax = std::min(xmax, e.xmax);
			ymin = std::max(ymin, e.ymin);
			ymax = std::min(ymax, e.ymax);
		}

		void unite(SpatExtent e) {
			xmin = std::min(xmin, e.xmin);
			xmax = std::max(xmax, e.xmax);
			ymin = std::min(ymin, e.ymin);
			ymax = std::max(ymax, e.ymax);
		}

		std::vector<double> asVector() {
			std::vector<double> e = { xmin, xmax, ymin, ymax };
			return(e);
		}

		std::vector<std::vector<double>> asPoints() {
			std::vector<std::vector<double>> pts(2, std::vector<double>(4));
			pts[0][0] = xmin;
			pts[1][0] = ymin;
			pts[0][1] = xmin;
			pts[1][1] = ymax;
			pts[0][2] = xmax;
			pts[1][2] = ymax;
			pts[0][3] = xmax;
			pts[1][3] = ymin;
			return(pts);
		}

		bool valid() {
			return ((xmax >= xmin) && (ymax >= ymin));
		}

		bool compare(SpatExtent e, std::string oper, double tolerance);

		SpatExtent round(int n);
		SpatExtent floor();
		SpatExtent ceil();
};



class SpatSRS {
	public:
		std::string proj4, wkt;
		bool set(std::string txt, std::string &msg);

/*
#ifdef useGDAL	
		bool set(OGRSpatialReference *poSRS, std::string &msg);
#endif		
*/
		std::string get(std::string x) {
			return (x == "proj4" ? proj4 : wkt); 
		}

		std::string get_prj() {
			return proj4;
		}

		bool is_equal(SpatSRS x) {
			return (proj4 == x.proj4);
		}

		bool is_empty() {
			return (wkt == "");
		}

		bool is_lonlat() {
			bool b1 = proj4.find("longlat") != std::string::npos;
			bool b2 = proj4.find("epsg:4326") != std::string::npos;
			return (b1 | b2);
		}

		bool could_be_lonlat(SpatExtent e) {
			bool b = is_lonlat();
			if ((!b) & is_empty()) {
				if ((e.xmin >= -180.1) & (e.xmax <= 180.1) & (e.ymin >= -90.1) & (e.ymax <= 90.1)) {
					b = true;
				}
			}
			return b;
		}

		bool is_global_lonlat(SpatExtent e) {
			if (could_be_lonlat(e)) {
                if (std::abs(e.xmax - e.xmin - 360) < 0.001) {
                    return true;
                }
				//double halfres = xres()/ 180;
				//if (((e.xmin - halfres) <= -180) && ((e.xmax + halfres) >= 180)) {
				//	return true;
				//}
			}
			return false;
		}
};

