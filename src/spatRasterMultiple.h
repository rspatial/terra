// Copyright (c) 2018-2023  Robert J. Hijmans
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

#include "spatRaster.h"


// A collection of (perhaps non matching) SpatRasters 
class SpatRasterCollection {
	public:
		virtual ~SpatRasterCollection(){}	
		SpatMessages msg;
		SpatRasterCollection deepCopy();
		void setError(std::string s);
		void addWarning(std::string s);
		bool has_error();
		bool has_warning();
		std::vector<std::string> getWarnings();
		std::string getError();
	
		std::vector<SpatRaster> ds;
//		SpatExtent extent;
		std::vector<std::string> names;
		SpatRasterCollection() {};
		SpatRasterCollection(std::string fname, std::vector<int> ids, bool useids);
//		void setExtent();
		SpatExtent getExtent();
		
		SpatRasterCollection(size_t n);
		size_t size();
		bool empty();
		void resize(size_t n);
		void push_back(SpatRaster r, std::string name);
		void erase(size_t i); 
		SpatRasterCollection crop(SpatExtent e, std::string snap, bool expand, std::vector<unsigned> use, SpatOptions &opt);
		SpatRasterCollection cropmask(SpatVector v, std::string snap, bool touches, bool expand, std::vector<unsigned> use, SpatOptions &opt);
		std::vector<int> getValueType(bool unique);

		SpatRaster merge(bool first, bool narm, SpatOptions &opt);
		SpatRaster morph(SpatRaster &x, SpatOptions &opt);
		SpatRaster mosaic(std::string fun, SpatOptions &opt);
		SpatRaster summary(std::string fun, SpatOptions &opt);
		std::vector<unsigned> dims();
		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> nms);
		std::vector<std::string> filenames();
		
};

// A class for "sub-datasets" 
class SpatRasterStack {
	public:
		virtual ~SpatRasterStack(){}
		SpatRasterStack deepCopy();
		SpatMessages msg;
		void setError(std::string s);
		void addWarning(std::string s);
		bool has_error();
		bool has_warning();
		std::vector<std::string> getWarnings();
		std::string getError();

		std::vector<SpatRaster> ds;
		std::vector<std::string> names;
		std::vector<std::string> long_names;
		std::vector<std::string> units;
		SpatRasterStack() {};
		SpatRasterStack(std::string fname, std::vector<int> ids, bool useids);
		SpatRasterStack(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn=false);
		std::vector<std::vector<std::vector<double>>> extractXY(std::vector<double> &x, std::vector<double> &y, std::string method);
		std::vector<std::vector<std::vector<double>>> extractCell(std::vector<double> &cell);
		std::vector<std::vector<std::vector<std::vector<double>>>> extractVector(SpatVector v, bool touches, std::string method, SpatOptions &opt);

		std::vector<double> resolution();
		SpatExtent getExtent();

		std::vector<std::string> get_names();
		void set_names(std::vector<std::string> nms);
		std::vector<std::string> get_longnames();
		void set_longnames(std::vector<std::string> nms);
		std::vector<std::string> get_units();
		void set_units(std::vector<std::string> u);
		std::vector<std::string> filenames();
		bool readStart();
		bool readStop();
		unsigned nsds();
		unsigned nrow();
		unsigned ncol();
		std::vector<unsigned> nlyr();

		std::string getSRS(std::string s);
		bool push_back(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn); 
		size_t size();
		bool empty();
		void resize(size_t n);
		void erase(size_t i);


		SpatRaster getsds(size_t i);
		SpatRasterStack subset(std::vector<unsigned> x);

		SpatRasterStack crop(SpatExtent e, std::string snap, bool expand, SpatOptions &opt);
		void replace(unsigned i, SpatRaster x);
		SpatRaster collapse();
		SpatRaster summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt);
		SpatRaster summary(std::string fun, bool narm, SpatOptions &opt);
};


