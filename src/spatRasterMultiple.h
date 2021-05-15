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

#include "spatRaster.h"


// A collection of (perhaps non matching) SpatRasters 
class SpatRasterCollection {
	public:
		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool has_error() { return msg.has_error; }
		bool has_warning() { return msg.has_warning; }
		std::string getWarnings() { return msg.getWarnings(); }
		std::string getError() { return msg.getError(); }
	
		std::vector<SpatRaster> ds;
		SpatRasterCollection() {};
		SpatRasterCollection(size_t n) { ds.resize(n); };
		size_t size() { return ds.size(); }
		void resize(size_t n) { ds.resize(n); }
		void push_back(SpatRaster r) { ds.push_back(r); };
		void erase(size_t i) { 
			if (i < ds.size()) {
				ds.erase(ds.begin()+i);
			}
		}
		SpatRaster merge(SpatOptions &opt);
		SpatRaster mosaic(std::string fun, SpatOptions &opt);
		SpatRaster summary(std::string fun, SpatOptions &opt);
		
};

// A class for "sub-datasets" 
class SpatRasterStack {
	public:
		SpatMessages msg;
		void setError(std::string s) { msg.setError(s); }
		void addWarning(std::string s) { msg.addWarning(s); }
		bool has_error() { return msg.has_error; }
		bool has_warning() { return msg.has_warning; }
		std::string getWarnings() { return msg.getWarnings();}
		std::string getError() { return msg.getError();}

		std::vector<SpatRaster> ds;
		std::vector<std::string> names;
		std::vector<std::string> long_names;
		std::vector<std::string> units;
		SpatRasterStack() {};
		SpatRasterStack(std::string fname, std::vector<int> ids, bool useids);
		SpatRasterStack(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn=false) { 
			push_back(r, name, longname, unit, warn); 
		};

		std::vector<std::vector<std::vector<double>>> extractXY(std::vector<double> &x, std::vector<double> &y, std::string method);
		std::vector<std::vector<std::vector<double>>> extractCell(std::vector<double> &cell);
		std::vector<std::vector<std::vector<std::vector<double>>>> extractVector(SpatVector v, bool touches, std::string method);

		std::vector<std::string> get_names() {
			return names;
		};
		void set_names(std::vector<std::string> nms) {
			if (nms.size() == ds.size()) {
				names = nms;
			}
		}
		std::vector<std::string> get_longnames() {
			return long_names;
		};
		void set_longnames(std::vector<std::string> nms) {
			if (nms.size() == ds.size()) {
				long_names = nms;
			}
		}
		std::vector<std::string> get_units() {
			return units;
		};
		void set_units(std::vector<std::string> u) {
			if (u.size() == ds.size()) {
				units = u;
			}
		}

		bool readStart() {
			for (auto& x : ds) { if (!x.readStart()) return false; }
			return true;
		}
			
		bool readStop() {
			for (auto& x : ds) { if (!x.readStop()) return false; }
			return true;
		}
	
		unsigned nsds() {
			return ds.size();
		}
		unsigned nrow() {
			if (ds.size() > 0) {
				return ds[0].nrow();
			} else {
				return 0;
			}
		}
		unsigned ncol() {
			if (ds.size() > 0) {
				return ds[0].ncol();
			} else {
				return 0;
			}
		}

		std::vector<unsigned> nlyr() {
			std::vector<unsigned> out;
			if (ds.size() > 0) {
				out.reserve(ds.size());
				for (size_t i=0; i<ds.size(); i++) {
					out.push_back(ds[i].nlyr());
				}
			} 
			return out;
		}

		std::string getSRS(std::string s) {
			if (ds.size() > 0) {
				return ds[0].getSRS(s);
			} else {
				return "";
			}
		}
		
		bool push_back(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn) { 
			if (ds.size() > 0) {
				if (!r.compare_geom(ds[0], false, false, true, true, true, false)) {
//				if (!ds[0].compare_geom(r, false, false, true, true, false, false)) {
					if (warn) {
						addWarning(r.msg.getError() +" (" + name + ")");
						return true;
					} else {
						setError(r.msg.getError() +" (" + name + ")");
						return false;
					}
				}
			}
			ds.push_back(r);
			names.push_back(name);
			long_names.push_back(longname);
			units.push_back(unit);
			return true;
		};
		
		size_t size() { return ds.size(); }
		void resize(size_t n) { 
			if (n < ds.size()) {
				ds.resize(n); 
				names.resize(n);
				long_names.resize(n);
				units.resize(n);
			}
		}
		void erase(size_t i) { 
			if (i < ds.size()) {
				ds.erase(ds.begin()+i); 
				names.erase(names.begin()+i);
				long_names.erase(long_names.begin()+i);
				units.erase(units.begin()+i);
			}	
		}


		SpatRaster getsds(size_t i) {
			if (i < ds.size()) {
				return(ds[i]); 
			} else {
				SpatRaster out;
				out.setError("invalid index");
				return out;
			}
		}
		SpatRasterStack subset(std::vector<unsigned> x) {
			SpatRasterStack out;
			for (size_t i=0; i<x.size(); i++) {
				if (x[i] < ds.size()) {
					out.push_back(ds[x[i]], names[i], long_names[i], units[i], true);
				} 				
			} 
			return out;
		}
		
		void replace(unsigned i, SpatRaster x) {
			if (i > (ds.size()-1)) {
				setError("invalid index");
				return;				
			}
			if (ds.size() == 0) {
				setError("cannot replace on empty stack");
				return;
			}
			if (!ds[0].compare_geom(x, false, false, true, true, false, false)) {
				setError("extent does not match");
				return;
			}
			
			ds[i] = x;
			names[i] = x.getNames()[0];
			long_names[i] = x.getLongSourceNames()[0];
			units[i] = x.getUnit()[0];
		}
		
		SpatRaster collapse() {
			SpatRaster out;

			if (ds.size() > 0) {
				out = ds[0];
				for (size_t i=1; i<ds.size(); i++) {
					for (size_t j=0; j<ds[i].source.size(); j++) {
						out.source.push_back(ds[i].source[j]);
					}
				}
			} 
			return out;
		}
		
		SpatRaster summary_numb(std::string fun, std::vector<double> add, bool narm, SpatOptions &opt);
		SpatRaster summary(std::string fun, bool narm, SpatOptions &opt);
};


