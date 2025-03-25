// Copyright (c) 2018-2025  Robert J. Hijmans
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
#include "string_utils.h"
#include "file_utils.h"


SpatRasterCollection SpatRasterCollection::deepCopy() { return *this; }

SpatRasterCollection::SpatRasterCollection(size_t n) { ds.resize(n); };

size_t SpatRasterCollection::size() { return ds.size(); }
bool SpatRasterCollection::empty() { return ds.empty(); }

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
	if (ds.empty()) {
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


std::string SpatRasterCollection::make_vrt(std::vector<std::string> options, bool reverse, SpatOptions &opt) {

	std::string outfile = opt.get_filename();
	if (outfile.empty()) {
		outfile = tempFile(opt.get_tempdir(), opt.tmpfile, ".vrt");
	} else if (file_exists(outfile) && (!opt.get_overwrite())) {
		setError("output file exists. You can use 'overwrite=TRUE' to overwrite it");
		return("");
	}
	opt.set_filenames({outfile});

	std::vector<std::string> ff;
	ff.reserve(size());
	
	SpatOptions xopt(opt);
	for (size_t i=0; i<size(); i++) {
//		if (!ds[i].hasValues()) continue;
		std::vector<std::string> f = ds[i].filenames();
		if ((ds[i].nsrc() == 1) && f[0] != "") {
			ff.push_back(f[0]);
		} else {
			std::string tmpf = tempFile(xopt.get_tempdir(), xopt.tmpfile, "_temp_raster.tif");
			xopt.set_filenames({tmpf});
			SpatRaster out = ds[i].writeRaster(xopt);
			if (out.hasError()) {
				setError(out.getError());
				return "";
			}
			ff.push_back(tmpf);
		}
	}
	SpatRaster tmp;
	if (reverse) std::reverse(ff.begin(), ff.end());
	return tmp.make_vrt(ff, options, opt);
}


void SpatRasterCollection::readBlock(SpatRaster &r, std::vector<std::vector<double>> &v, BlockSize bs, size_t i, std::vector<size_t> use, SpatOptions opt){

	if ((bs.row[i] + bs.nrows[i]) > r.nrow()) {
		setError("invalid rows/columns");
		return;
	}
	if (bs.nrows[i]==0) {
		return;
	}
	SpatExtent re = r.getExtent();
	double yres = r.yres();
	double ymx = re.ymax - bs.row[i] * yres;
	double ymn = re.ymax - (bs.row[i] + bs.nrows[i]) * yres;
	SpatExtent e = {re.xmin, re.xmax, ymn, ymx};
	SpatRasterCollection x = crop(e, "near", true, use, opt);
	if (x.hasError()) {
		setError(x.getError());
		return;
	}
	v.resize(x.size());
	for (size_t i=0; i< x.size(); i++) {
		x.ds[i].readValues(v[i], 0, x.ds[i].nrow(), 0, x.ds[i].ncol());
	}
}



SpatRasterCollection SpatRasterCollection::crop(SpatExtent e, std::string snap, bool expand, std::vector<size_t> use, SpatOptions &opt) {

	SpatRasterCollection out;
	if ( !e.valid() ) {
		out.setError("invalid extent");
		return out;
	}
	if (!e.valid_notempty()) {
		out.setError("cannot crop with an empty extent");
		return out;
	}
	SpatOptions ops(opt);
	if (use.empty()) {
		for (size_t i=0; i<size(); i++) {
			SpatExtent xe = e.intersect(ds[i].getExtent());
			if (xe.valid_notempty()) {
				SpatRaster r = ds[i].crop(e, snap, expand, ops);
				if (r.hasError()) {
					out.setError(r.getError());
					return out;
				} 
				out.push_back(r, "");
			}
		}
	} else {
		for (size_t i=0; i<use.size(); i++) {
			SpatExtent xe = e.intersect(ds[use[i]].getExtent());
			if (xe.valid_notempty()) {
				SpatRaster r = ds[use[i]].crop(e, snap, expand, ops);
				if (r.hasError()) {
					out.setError(r.getError());
					return out;
				}
				out.push_back(r, "");
			}
		}
	}
	return out;
}


SpatRasterCollection SpatRasterCollection::cropmask(SpatVector v, std::string snap, bool touches, bool expand, std::vector<size_t> use, SpatOptions &opt) {
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
	if (use.empty()) {
		for (size_t i=0; i<size(); i++) {
			SpatExtent xe = e.intersect(ds[i].getExtent());
			if (xe.valid()) {
				SpatRaster x = ds[i].cropmask(v, snap, touches, expand, ops);
				out.push_back(x.source[0], names[i]);
			}
		}
	} else {
		for (size_t i=0; i<use.size(); i++) {
			SpatExtent xe = e.intersect(ds[use[i]].getExtent());
			if (xe.valid()) {
				SpatRaster r = ds[use[i]].cropmask(v, snap, touches, expand, ops);
				out.push_back(r.source[0], names[use[i]]);
			}
		}
	}
	return out;
}

std::vector<size_t> SpatRasterCollection::dims() {
	size_t n = ds.size();
	size_t n2 = 2 * n;
	std::vector<size_t> out(n * 3);
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
		


bool SpatRasterCollection::addTag(std::string name, std::string value, std::string domain) {
	lrtrim(name);
	lrtrim(value);
	if (value == "") {
		return removeTag(name, domain);
	} else if (name != "") {
		if (tags.size() == 0) {
			tags.resize(3);
		}
		tags[0].push_back(domain);
		tags[1].push_back(name);
		tags[2].push_back(value);
		std::sort(tags.begin(), tags.end());
		tags.erase(std::unique(tags.begin(), tags.end()), tags.end());
		return true;
	} 
	return false;
}


bool SpatRasterCollection::removeTag(std::string name, std::string domain) {
	if (tags.size() == 0) return true;
	for (size_t i =0; i<tags[0].size(); i++) {
		if ((tags[0][i] == domain) && (tags[1][i] == name)) {
			tags[0].erase(tags[0].begin()+i);
			tags[1].erase(tags[1].begin()+i);
			tags[2].erase(tags[2].begin()+i);
			return true;
		}
	}
	return false;
}


std::string SpatRasterCollection::getTag(std::string name, std::string domain) {
	for (size_t i =0; i<tags[0].size(); i++) {
		if ((tags[0][i] == domain) && (tags[1][i] == name)) {
			return tags[2][i];
		}
	}
	return "";
}

std::vector<std::vector<std::string>> SpatRasterCollection::getTags() {
	return tags;
}



/////////////////////////////////////////////////

SpatRasterStack SpatRasterStack::deepCopy() { return *this; }


SpatRasterStack::SpatRasterStack(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn) { 
	push_back(r, name, longname, unit, warn); 
};


std::vector<double> SpatRasterStack::resolution() {
	if (ds.empty()) {
		return {NAN, NAN};
	} else {
		return ds[0].resolution();
	}
}
SpatExtent SpatRasterStack::getExtent() {
	if (ds.empty()) {
		return SpatExtent();
	} else {
		return ds[0].getExtent();
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

void SpatRasterStack::set_layernames(std::vector<std::string> nms, long id) {
	if (id < 0) {
		for (size_t i=0; i<ds.size(); i++) {
			if (ds[i].nlyr() != nms.size()) {
				setError("length of names does not match the number of layers");
			} else {
				ds[i].setNames(nms);
			}
		}
	} else {
		if (ds[id].nlyr() != nms.size()) {
			setError("length of names does not match the number of layers");
		} else {
			ds[id].setNames(nms);
		}
	}
}

std::vector<std::vector<std::string>> SpatRasterStack::get_layernames() {
	size_t nd = ds.size();
	std::vector<std::vector<std::string>> out(nd);
	for (size_t i = 0; i<nd; i++) {
		out[i] = ds[i].getNames();
	}		
	return out;
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
	
bool SpatRasterStack::readAll() {
  for (auto& x : ds) { if (!x.readAll()) return false; }
  return true;
}

size_t SpatRasterStack::nsds() {
	return ds.size();
}
size_t SpatRasterStack::nrow() {
	if (ds.empty()) {
		return 0;
	} else {
		return ds[0].nrow();
	}
}
size_t SpatRasterStack::ncol() {
	if (ds.empty()) {
		return 0;
	} else {
		return ds[0].ncol();
	}
}
std::vector<size_t> SpatRasterStack::nlyr() {
	std::vector<size_t> out;
	if (!ds.empty()) {
		out.reserve(ds.size());
		for (size_t i=0; i<ds.size(); i++) {
			out.push_back(ds[i].nlyr());
		}
	} 
	return out;
}

std::string SpatRasterStack::getSRS(std::string s) {
	if (ds.empty()) {
		return "";
	} else {
		return ds[0].getSRS(s);
	}
}
		
bool SpatRasterStack::push_back(SpatRaster r, std::string name, std::string longname, std::string unit, bool warn) { 
	if (!ds.empty()) {
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
bool SpatRasterStack::empty() { return ds.empty(); }
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

SpatRasterStack SpatRasterStack::subset(std::vector<size_t> x) {
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
		if (hasError()) {
			return(out);
		}
	}
	return out;
}

void SpatRasterStack::replace(size_t i, SpatRaster x, bool setname) {
	if (i > (ds.size()-1)) {
		setError("invalid index");
		return;				
	}
	if (ds.empty()) {
		setError("cannot replace on empty stack");
		return;
	}
	if (!ds[0].compare_geom(x, false, false, true, true, false, false)) {
		setError("extent does not match");
		return;
	}
	
	ds[i] = x;
//  for clause for #1604
	if (setname) {
		names[i] = x.getNames()[0];
		long_names[i] = x.getLongSourceNames()[0];
		units[i] = x.getUnit()[0];
	}
}

SpatRaster SpatRasterStack::collapse() {
	SpatRaster out;

	if (!ds.empty()) {
		out = ds[0];
		for (size_t i=1; i<ds.size(); i++) {
			for (size_t j=0; j<ds[i].source.size(); j++) {
				out.source.push_back(ds[i].source[j]);
			}
		}
	} 
	out.user_tags = tags;
	return out;
}



bool SpatRasterStack::addTag(std::string name, std::string value, std::string domain) {
	lrtrim(name);
	lrtrim(value);
	if (value == "") {
		return removeTag(name, domain);
	} else if (name != "") {
		if (tags.size() == 0) {
			tags.resize(3);
		}
		tags[0].push_back(domain);
		tags[1].push_back(name);
		tags[2].push_back(value);
		std::sort(tags.begin(), tags.end());
		tags.erase(std::unique(tags.begin(), tags.end()), tags.end());
		return true;
	} 
	return false;
}


bool SpatRasterStack::removeTag(std::string name, std::string domain) {
	if (tags.size() == 0) return true;
	for (size_t i =0; i<tags[0].size(); i++) {
		if ((tags[0][i] == domain) && (tags[1][i] == name)) {
			tags[0].erase(tags[0].begin()+i);
			tags[1].erase(tags[1].begin()+i);
			tags[2].erase(tags[2].begin()+i);
			return true;
		}
	}
	return false;
}


std::string SpatRasterStack::getTag(std::string name, std::string domain) {
	for (size_t i =0; i<tags[0].size(); i++) {
		if ((tags[0][i] == domain) && (tags[1][i] == name)) {
			return tags[2][i];
		}
	}
	return "";
}

std::vector<std::vector<std::string>> SpatRasterStack::getTags() {
	return tags;
}

