// Copyright (c) 2018-2022  Robert J. Hijmans
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

#include "spatRasterMultiple.h"

SpatRasterCollection SpatRasterCollection::deepCopy() { return *this; }
void SpatRasterCollection::setError(std::string s) { msg.setError(s); }
void SpatRasterCollection::addWarning(std::string s) { msg.addWarning(s); }
bool SpatRasterCollection::has_error() { return msg.has_error; }
bool SpatRasterCollection::has_warning() { return msg.has_warning; }
std::vector<std::string> SpatRasterCollection::getWarnings() { return msg.getWarnings(); }
std::string SpatRasterCollection::getError() { return msg.getError(); }

SpatRasterCollection::SpatRasterCollection(size_t n) { ds.resize(n); };

size_t SpatRasterCollection::size() { return ds.size(); }

void SpatRasterCollection::resize(size_t n) { ds.resize(n); }

void SpatRasterCollection::push_back(SpatRaster r, std::string name) { 

/*
	if (ds.size() == 0) {
		extent = r.getExtent();
	} else {
		extent.unite(r.getExtent());
	}
*/
	ds.push_back(r);
	names.push_back(name);
	
}

SpatExtent SpatRasterCollection::getExtent() { 
	SpatExtent e;
	if (ds.size() == 0) {
		e = SpatExtent();
	} else {
		e = ds[0].getExtent();
		for (size_t i=1; i<ds.size(); i++) {
			e.unite(ds[i].getExtent());
		}
	}
	return e;
}


std::vector<int> SpatRasterCollection::getValueType(bool unique) {
	std::vector<int> d;
	
	for (size_t i=0; i<ds.size(); i++) {
		std::vector<int> dd = ds[i].getValueType(false);		
		d.insert(d.end(), dd.begin(), dd.end());
	}
	if (unique) {
		std::sort(d.begin(), d.end());
		d.erase(std::unique(d.begin(), d.end()), d.end());	
	}
	return(d);
}



/*
void SpatRasterCollection::setExtent() { 
	if (ds.size() == 0) {
		extent = SpatExtent();
		return;
	} else {
		extent = ds[0].getExtent();
	}
	for (size_t i=1; i<ds.size(); i++) {
		extent.unite(ds[i].getExtent());
	}
}

SpatExtent SpatRasterCollection::getExtent() { 
	return extent;
}
*/

void SpatRasterCollection::erase(size_t i) { 
	if (i < ds.size()) {
		ds.erase(ds.begin()+i);
//		setExtent();
	}
}

SpatRasterCollection SpatRasterCollection::crop(SpatExtent e, std::string snap, bool expand, std::vector<unsigned> use, SpatOptions &opt) {


	SpatRasterCollection out;
	if ( !e.valid() ) {
		out.setError("invalid extent");
		return out;
	}
	if ((e.xmin == e.xmax) && (e.ymin == e.ymax)) {
		out.setError("cannot crop with an empty extent");
		return out;
	}
	SpatOptions ops(opt);
	if (use.size() > 0) {
		for (size_t i=0; i<use.size(); i++) {
			SpatExtent xe = e.intersect(ds[use[i]].getExtent());
			if (xe.valid()) {
				SpatRaster r = ds[use[i]];
				r = r.crop(e, snap, expand, ops);
				out.push_back(r, "");
			}
		}
	} else {
		for (size_t i=0; i<size(); i++) {
			SpatExtent xe = e.intersect(ds[i].getExtent());
			if (xe.valid()) {
				out.push_back(ds[i].crop(e, snap, expand, ops), "");
			}
		}
	}
	return out;
}


SpatRasterCollection SpatRasterCollection::cropmask(SpatVector v, std::string snap, bool touches, bool expand, std::vector<unsigned> use, SpatOptions &opt) {
	SpatRasterCollection out;

	SpatExtent e = v.extent;
	if ( !e.valid() ) {
		out.setError("invalid extent");
		return out;
	}
	if ((e.xmin == e.xmax) && (e.ymin == e.ymax)) {
		out.setError("cannot crop with an empty extent");
		return out;
	}
	SpatOptions ops(opt);
	if (use.size() > 0) {
		for (size_t i=0; i<use.size(); i++) {
			SpatExtent xe = e.intersect(ds[use[i]].getExtent());
			if (xe.valid()) {
				SpatRaster r = ds[use[i]].cropmask(v, snap, touches, expand, ops);
				out.push_back(r.source[0], names[use[i]]);
			}
		}
	} else {
		for (size_t i=0; i<size(); i++) {
			SpatExtent xe = e.intersect(ds[i].getExtent());
			if (xe.valid()) {
				SpatRaster x = ds[i].cropmask(v, snap, touches, expand, ops);
				out.push_back(x.source[0], names[i]);
			}
		}
	}
	return out;
}

std::vector<unsigned> SpatRasterCollection::dims() {
	size_t n = ds.size();
	size_t n2 = 2 * n;
	std::vector<unsigned> out(n * 3);
	for (size_t i=0; i<n; i++) {
		out[i]    = ds[i].nrow();
		out[i+n]  = ds[i].ncol();
		out[i+n2] = ds[i].nlyr();
	}
	return out;
};

std::vector<std::string> SpatRasterCollection::get_names() {
	return names;
};
void SpatRasterCollection::set_names(std::vector<std::string> nms) {
	if (nms.size() == ds.size()) {
		names = nms;
	}
}

std::vector<std::string> SpatRasterCollection::filenames() {
	size_t n =0;
	for (size_t i=0; i<ds.size(); i++) { 
		n += ds[i].nlyr();
	}
	std::vector<std::string> names;
	names.reserve(n);
	for (size_t i=0; i<ds.size(); i++) { 
		std::vector<std::string> n = ds[i].filenames();
		names.insert(names.end(), n.begin(), n.end());
	}
	return names;
};
		


SpatRasterStack SpatRasterStack::deepCopy() { return *this; }
void SpatRasterStack::setError(std::string s) { msg.setError(s); }
void SpatRasterStack::addWarning(std::string s) { msg.addWarning(s); }
bool SpatRasterStack::has_error() { return msg.has_error; }
bool SpatRasterStack::has_warning() { return msg.has_warning; }
std::vector<std::string> SpatRasterStack::getWarnings() { return msg.getWarnings();}
std::string SpatRasterStack::getError() { return msg.getError();}

SpatRasterStack::SpatRasterStack(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn) { 
	push_back(r, name, longname, unit, warn); 
};


std::vector<double> SpatRasterStack::resolution() {
	if (ds.size() > 0) {
		return ds[0].resolution();
	} else {
		return {NAN, NAN};
	}
}
SpatExtent SpatRasterStack::getExtent() {
	if (ds.size() > 0) {
		return ds[0].getExtent();
	} else {
		return SpatExtent();
	}
}

std::vector<std::string> SpatRasterStack::get_names() {
	return names;
};
void SpatRasterStack::set_names(std::vector<std::string> nms) {
	if (nms.size() == ds.size()) {
		names = nms;
	}
}
std::vector<std::string> SpatRasterStack::get_longnames() {
	return long_names;
};
void SpatRasterStack::set_longnames(std::vector<std::string> nms) {
	if (nms.size() == ds.size()) {
		long_names = nms;
	}
}
std::vector<std::string> SpatRasterStack::get_units() {
	return units;
};
void SpatRasterStack::set_units(std::vector<std::string> u) {
	if (u.size() == ds.size()) {
		units = u;
	}
}
std::vector<std::string> SpatRasterStack::filenames() {
	size_t n =0;
	for (size_t i=0; i<ds.size(); i++) { 
		n += ds[i].nlyr();
	}
	std::vector<std::string> names;
	names.reserve(n);
	for (size_t i=0; i<ds.size(); i++) { 
		std::vector<std::string> n = ds[i].filenames();
		names.insert(names.end(), n.begin(), n.end());
	}
	return names;
};

bool SpatRasterStack::readStart() {
	for (auto& x : ds) { if (!x.readStart()) return false; }
	return true;
}
			
bool SpatRasterStack::readStop() {
	for (auto& x : ds) { if (!x.readStop()) return false; }
	return true;
}
	
unsigned SpatRasterStack::nsds() {
	return ds.size();
}
unsigned SpatRasterStack::nrow() {
	if (ds.size() > 0) {
		return ds[0].nrow();
	} else {
		return 0;
	}
}
unsigned SpatRasterStack::ncol() {
	if (ds.size() > 0) {
		return ds[0].ncol();
	} else {
		return 0;
	}
}
std::vector<unsigned> SpatRasterStack::nlyr() {
	std::vector<unsigned> out;
	if (ds.size() > 0) {
		out.reserve(ds.size());
		for (size_t i=0; i<ds.size(); i++) {
			out.push_back(ds[i].nlyr());
		}
	} 
	return out;
}

std::string SpatRasterStack::getSRS(std::string s) {
	if (ds.size() > 0) {
		return ds[0].getSRS(s);
	} else {
		return "";
	}
}
		
bool SpatRasterStack::push_back(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn) { 
	if (ds.size() > 0) {
		if (!r.compare_geom(ds[0], false, false, true, true, true, false)) {
//		if (!ds[0].compare_geom(r, false, false, true, true, false, false)) {
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
	if (r.hasWarning()) {
		for (size_t i=0; i<r.msg.warnings.size(); i++) {
			addWarning(r.msg.warnings[i]);
		}
	}
	if (r.hasError()) {
		setError(r.msg.getError());
		return false;
	}
	return true;
};
		
size_t SpatRasterStack::size() { return ds.size(); }
void SpatRasterStack::resize(size_t n) { 
	if (n < ds.size()) {
		ds.resize(n); 
		names.resize(n);
		long_names.resize(n);
		units.resize(n);
	}
}
void SpatRasterStack::erase(size_t i) { 
	if (i < ds.size()) {
		ds.erase(ds.begin()+i); 
		names.erase(names.begin()+i);
		long_names.erase(long_names.begin()+i);
		units.erase(units.begin()+i);
	}	
}


SpatRaster SpatRasterStack::getsds(size_t i) {
	if (i < ds.size()) {
		return(ds[i]); 
	} else {
		SpatRaster out;
		out.setError("invalid index");
		return out;
	}
}

SpatRasterStack SpatRasterStack::subset(std::vector<unsigned> x) {
	SpatRasterStack out;
	for (size_t i=0; i<x.size(); i++) {
		size_t j = x[i];
		if (j < ds.size()) {
			out.push_back(ds[j], names[j], long_names[j], units[j], true);
		} 				
	} 
	return out;
}

SpatRasterStack SpatRasterStack::crop(SpatExtent e, std::string snap, bool expand, SpatOptions &opt) {
	SpatRasterStack out;
	std::vector<std::string> ff = opt.get_filenames();
	if (ff.size() != ds.size()) {
		opt.set_filenames({""});
		opt.ncopies *= ds.size();
	}
	for (size_t i=0; i<ds.size(); i++) {
		out.push_back(ds[i].crop(e, snap, expand, opt), names[i], long_names[i], units[i], true);
		if (has_error()) {
			return(out);
		}
	}
	return out;
}

void SpatRasterStack::replace(unsigned i, SpatRaster x) {
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

SpatRaster SpatRasterStack::collapse() {
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

