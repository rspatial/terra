// Copyright (c) 2018  Robert J. Hijmans
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

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


class SpatMessages {
	public:
		bool success = true;
		bool has_error = false;
		bool has_warning = false;
		std::string error;
		std::vector<std::string> warnings;
		
		void setError(std::string s) {
			has_error = true;
			error = s;
			success = false;
		}
		void addWarning(std::string s) {
			has_warning = true;
			warnings.push_back(s);
		}
};



class SpatOptions {
	private:
		std::string temp_dir = "";
		std::string file_name = "";
		std::string data_type = "FLT4S";
		std::string file_type = "GTiff"; 
		bool over_write = false;
		double mem_frac = 0.6;
		std::vector<std::string> gdal_options;
	public:
		bool overwrite = false;
		SpatOptions();
		SpatOptions(const SpatOptions &opt);
		SpatOptions deepcopy(const SpatOptions &opt);
		std::string get_filename();
		std::string get_datatype();
		std::string get_filetype();
		std::string get_tempdir();
		bool get_overwrite();
		void set_filename(std::string d);
		void set_datatype(std::string d);
		void set_filetype(std::string d);
		void set_tempdir(std::string d);
		void set_overwrite(bool d);

		void set_memfrac(double d);
		double get_memfrac();
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
		
		void intersect(SpatExtent e) { 
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
			std::vector<double> e(4);
			e[0] = xmin; e[1] = xmax; e[2] = ymin; e[3] = ymax; 
			return(e);
		}
			
		bool valid() {
			return ((xmax > xmin) && (ymax > ymin));
		}
};





