// Copyright (c) 2018-2021  Robert J. Hijmans
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

#include "spatDataframe.h"
#include <vector>
#include <algorithm>
#include <string>
#include "NA.h"
#include "string_utils.h"


SpatDataFrame::SpatDataFrame() {}

SpatDataFrame SpatDataFrame::skeleton() {
	SpatDataFrame out;
	out.names  = names;
	out.itype  = itype;
	out.iplace = iplace;
	out.dv = std::vector<std::vector<double>>(dv.size());
	out.iv = std::vector<std::vector<long>>(iv.size());
	out.sv = std::vector<std::vector<std::string>>(sv.size());
	out.bv = std::vector<std::vector<int8_t>>(bv.size());
	out.tv = std::vector<std::vector<SpatTime_t>>(tv.size());
	return out;
}


std::vector<double> SpatDataFrame::getD(unsigned i) {
	unsigned j = iplace[i];
	return dv[j];
}

double SpatDataFrame::getDvalue(unsigned i, unsigned j) {
	j = iplace[j];
	return dv[j][i];
}

std::vector<long> SpatDataFrame::getI(unsigned i) {
	unsigned j = iplace[i];
	return iv[j];
}


long SpatDataFrame::getIvalue(unsigned i, unsigned j) {
	j = iplace[j];
	return iv[j][i];
}


std::vector<std::string> SpatDataFrame::getS(unsigned i) {
	unsigned j = iplace[i];
	return sv[j];
}

std::string SpatDataFrame::getSvalue(unsigned i, unsigned j) {
	j = iplace[j];
	return sv[j][i];
}

std::vector<int8_t> SpatDataFrame::getB(unsigned i) {
	unsigned j = iplace[i];
	return bv[j];
}

int8_t SpatDataFrame::getBvalue(unsigned i, unsigned j) {
	j = iplace[j];
	return bv[j][i];
}

std::vector<SpatTime_t> SpatDataFrame::getT(unsigned i) {
	unsigned j = iplace[i];
	return tv[j];
}

SpatTime_t SpatDataFrame::getTvalue(unsigned i, unsigned j) {
	j = iplace[j];
	return tv[j][i];
}


SpatDataFrame SpatDataFrame::subset_rows(unsigned i) {
	std::vector<unsigned> r = { i }; 
	return subset_rows(r);
}


SpatDataFrame SpatDataFrame::subset_rows(std::vector<unsigned> range) { 

	SpatDataFrame out;
	unsigned nr = nrow();
	if (range.size() > 0) {
		for (int i = range.size(); i>0; i--) {
			if (range[i-1] > nr) {
				range.erase(range.begin() + i-1);
			}
		}
	}
	out.names = names;
	out.itype = itype;
	out.iplace = iplace;

	out.dv.resize(dv.size());
	out.iv.resize(iv.size());
	out.sv.resize(sv.size());
	out.bv.resize(bv.size());
	out.tv.resize(tv.size());
	out.reserve(range.size());
	
	for (size_t i=0; i < range.size(); i++) {
		for (size_t j=0; j < dv.size(); j++) {
			out.dv[j].push_back(dv[j][range[i]]);
		}
		for (size_t j=0; j < iv.size(); j++) {
			out.iv[j].push_back(iv[j][range[i]]);
		}
		for (size_t j=0; j < sv.size(); j++) {
			out.sv[j].push_back(sv[j][range[i]]);
		}
		for (size_t j=0; j < bv.size(); j++) {
			out.bv[j].push_back(bv[j][range[i]]);
		}
		for (size_t j=0; j < tv.size(); j++) {
			out.tv[j].push_back(tv[j][range[i]]);
		}
	}
	return out;
}


SpatDataFrame SpatDataFrame::subset_rows(std::vector<long> range) { 
	std::vector<unsigned> r(range.begin(), range.end());
	return subset_rows(r);
}


SpatDataFrame SpatDataFrame::subset_cols(unsigned i) {
	std::vector<unsigned> c = { i }; 
	return subset_cols(c);
}

SpatDataFrame SpatDataFrame::subset_cols(std::vector<unsigned> range) { 
	SpatDataFrame out;
	unsigned dcnt=0;
	unsigned icnt=0;
	unsigned scnt=0;
	unsigned bcnt=0;
	unsigned tcnt=0;
	for (size_t i=0; i < range.size(); i++) {
		if (range[i] < 0 || range[i] >= ncol()) {
			out.setError("invalid column");
			return out;
		}
		unsigned j = range[i];
		unsigned p = iplace[j];
		out.names.push_back(names[j]);
		if (itype[j] == 0) {
			out.dv.push_back(dv[p]);
			out.iplace.push_back(dcnt);
			out.itype.push_back(0);
			dcnt++;
		} else if (itype[j] == 1) {
			out.iv.push_back(iv[p]);
			out.iplace.push_back(icnt);
			out.itype.push_back(1);
			icnt++;
		} else if (itype[j] == 2) {
			out.sv.push_back(sv[p]);
			out.iplace.push_back(scnt);
			out.itype.push_back(2);
			scnt++;
		} else if (itype[j] == 3) {
			out.bv.push_back(bv[p]);
			out.iplace.push_back(bcnt);
			out.itype.push_back(3);
			bcnt++;
		} else {
			out.tv.push_back(tv[p]);
			out.iplace.push_back(tcnt);
			out.itype.push_back(4);
			tcnt++;
		}
	}
	return out;
}


unsigned SpatDataFrame::ncol() {
	return itype.size();
}


unsigned SpatDataFrame::nrow() {
	unsigned n;
	if (itype.size() == 0) {
		n = 0;
	} else {
		if (itype[0] == 0) {
			n = dv[0].size();
		} else if (itype[0] == 1) {
			n = iv[0].size();
		} else if (itype[0] == 2) {
			n = sv[0].size();
		} else if (itype[0] == 3) {
			n = bv[0].size();
		} else {
			n = tv[0].size();
		}
	}
	return n;
}




void SpatDataFrame::add_row() {
	for (size_t i=0; i < dv.size(); i++) {
		dv[i].push_back(NAN);
	}
	long longNA = NA<long>::value;
	for (size_t i=0; i < iv.size(); i++) {
		iv[i].push_back(longNA);
	}
	for (size_t i=0; i < sv.size(); i++) {
		sv[i].push_back(NAS);
	}
	for (size_t i=0; i < bv.size(); i++) {
		bv[i].push_back(2);
	}
	SpatTime_t timeNA = NA<SpatTime_t>::value;
	for (size_t i=0; i < tv.size(); i++) {
		tv[i].push_back(timeNA);
	}
}


void SpatDataFrame::add_rows(size_t n) {
	size_t s = nrow() + n;
	for (size_t i=0; i < dv.size(); i++) {
		dv[i].resize(s, NAN);
	}
	long longNA = NA<long>::value;
	for (size_t i=0; i < iv.size(); i++) {
		iv[i].resize(s, longNA);
	}
	for (size_t i=0; i < sv.size(); i++) {
		sv[i].resize(s, NAS);
	}
	for (size_t i=0; i < bv.size(); i++) {
		bv[i].resize(s, 2);
	}
	SpatTime_t timeNA = NA<SpatTime_t>::value;
	for (size_t i=0; i < tv.size(); i++) {
		tv[i].resize(s, timeNA);
	}
}


void SpatDataFrame::reserve(unsigned n) {
	for (size_t i=0; i<dv.size(); i++) {
		dv[i].reserve(n);
	}
	for (size_t i=0; i<iv.size(); i++) {
		iv[i].reserve(n);
	}
	for (size_t i=0; i<sv.size(); i++) {
		sv[i].reserve(n);
	}
	for (size_t i=0; i<bv.size(); i++) {
		bv[i].reserve(n);
	}
	for (size_t i=0; i<tv.size(); i++) {
		tv[i].reserve(n);
	}
}

void SpatDataFrame::resize_rows(unsigned n) {
	for (size_t i=0; i<dv.size(); i++) {
		dv[i].resize(n, NAN);
	}
	long longNA = NA<long>::value;
	for (size_t i=0; i<iv.size(); i++) {
		iv[i].resize(n, longNA);
	}
	for (size_t i=0; i<sv.size(); i++) {
		sv[i].resize(n, NAS);
	}
	for (size_t i=0; i<bv.size(); i++) {
		bv[i].resize(n, 2);
	}
	SpatTime_t timeNA = NA<SpatTime_t>::value;
	for (size_t i=0; i<tv.size(); i++) {
		tv[i].resize(n, timeNA);
	}
}

void SpatDataFrame::remove_rows(std::vector<unsigned> r) {
	if (r.size() == 0) return;
	//sort(r.begin(), r.end(), std::greater<unsigned>());
	sort(r.begin(), r.end());
	r.erase(std::unique(r.begin(), r.end()), r.end());
	std::reverse(r.begin(), r.end());

	for (size_t j=0; j<r.size(); j++) {
		for (size_t i=0; i<dv.size(); i++) {
			dv[i].erase(dv[i].begin() + r[j]);
		}
		for (size_t i=0; i<iv.size(); i++) {
			iv[i].erase(iv[i].begin() +r[j]);
		}
		for (size_t i=0; i<sv.size(); i++) {
			sv[i].erase(sv[i].begin() +r[j]);
		}
		for (size_t i=0; i<bv.size(); i++) {
			bv[i].erase(bv[i].begin() +r[j]);
		}
		for (size_t i=0; i<tv.size(); i++) {
			tv[i].erase(tv[i].begin() +r[j]);
		}
	}
}



void SpatDataFrame::resize_cols(unsigned n) {
	if (n < ncol()) {
		itype.resize(n);
		iplace.resize(n);
	} else {
		setError("you can only resize to fewer columns");
	}
}

// use template instead
bool SpatDataFrame::add_column(std::vector<double> x, std::string name) {
	unsigned nr = nrow();
	if ((nr != 0) & (nr != x.size())) return false; 
	iplace.push_back(dv.size());
	itype.push_back(0);
	names.push_back(name);
	dv.push_back(x);
	return true;
}


bool SpatDataFrame::remove_column(int i) {

	if ((i < 0) | ((size_t)i > ncol())) {
		return false;
	}
	size_t dtype = itype[i];
	size_t place = iplace[i];

	size_t ii = i;
	if (ii < (iplace.size()-1)) {
		for (size_t j=i+1; j<iplace.size(); j++) {
			if (itype[j] == dtype) {
				iplace[j]--;
			}
		}
	}

	names.erase(names.begin()+i);
	itype.erase(itype.begin()+i);
	iplace.erase(iplace.begin()+i);
	if (dtype == 0) {
		dv.erase(dv.begin()+place);
	} else if (dtype == 1) {
		iv.erase(iv.begin()+place);
	} else if (dtype == 2) {
		sv.erase(sv.begin()+place);
	} else if (dtype == 3) {
		bv.erase(bv.begin()+place);
	} else {
		tv.erase(tv.begin()+place);
	}
	return true;
}

bool SpatDataFrame::remove_column(std::string field) {
	int i = where_in_vector(field, names, false);
	return remove_column(i);
}



bool SpatDataFrame::add_column(std::vector<long> x, std::string name) {
	unsigned nr = nrow();
	if ((nr != 0) & (nr != x.size())) return false; 
	iplace.push_back(iv.size());
	itype.push_back(1);
	names.push_back(name);
	iv.push_back(x);
	return true;
}

bool SpatDataFrame::add_column(std::vector<int> x, std::string name) {
	std::vector<long> v(x.begin(), x.end());
	return add_column(v, name);
}


bool SpatDataFrame::add_column(std::vector<std::string> x, std::string name) {
	unsigned nr = nrow();
	if ((nr != 0) & (nr != x.size())) return false; 
	iplace.push_back(sv.size());
	itype.push_back(2);
	names.push_back(name);
	sv.push_back(x);
	return true;
}

bool SpatDataFrame::add_column(std::vector<int8_t> x, std::string name) {
	unsigned nr = nrow();
	if ((nr != 0) & (nr != x.size())) return false; 
	iplace.push_back(bv.size());
	itype.push_back(3);
	names.push_back(name);
	bv.push_back(x);
	return true;
}

bool SpatDataFrame::add_column_bool(std::vector<int> x, std::string name) {
	unsigned nr = nrow();
	if ((nr != 0) & (nr != x.size())) return false; 
	iplace.push_back(bv.size());
	itype.push_back(3);
	names.push_back(name);
	std::vector<int8_t> b;
	b.reserve(x.size());
	for (size_t i=0; i<x.size(); i++){
		if ((x[i] == 0)|| (x[i]==1)) {
			b.push_back(x[i]);
		} else {
			b.push_back(2);
		}
	}
	bv.push_back(b);
	return true;
}


bool SpatDataFrame::add_column(std::vector<SpatTime_t> x, std::string name) {
	unsigned nr = nrow();
	if ((nr != 0) & (nr != x.size())) return false; 
	iplace.push_back(tv.size());
	itype.push_back(4);
	names.push_back(name);
	tv.push_back(x);
	return true;
}


void SpatDataFrame::add_column(unsigned dtype, std::string name) {
	unsigned nr = nrow();
	if (dtype == 0) {
		std::vector<double> dins(nr, NAN);
		iplace.push_back(dv.size());
		dv.push_back(dins);
	} else if (dtype == 1) {
		long longNA = NA<long>::value;
		std::vector<long> iins(nr, longNA);
		iplace.push_back(iv.size());
		iv.push_back(iins);
	} else if (dtype == 2) {
		std::vector<std::string> sins(nr, NAS);
		iplace.push_back(sv.size());
		sv.push_back(sins);
	} else if (dtype == 3) {
		std::vector<int8_t> bins(nr, 2);
		iplace.push_back(bv.size());
		bv.push_back(bins);
	} else {
		SpatTime_t timeNA = NA<SpatTime_t>::value;
		std::vector<SpatTime_t> sins(nr, timeNA);
		iplace.push_back(tv.size());
		tv.push_back(sins);
	}
	itype.push_back(dtype);
	names.push_back(name);
}

bool SpatDataFrame::cbind(SpatDataFrame &x) {
	unsigned nc = x.ncol();
	std::vector<std::string> nms = x.names;
	for (size_t i=0; i<nc; i++) {
		if (x.itype[i] == 0) {
			std::vector<double> d = x.getD(i);
			if (!add_column(d, nms[i])) return false;
		} else if (x.itype[i] == 1) {
			std::vector<long> d = x.getI(i);
			if (!add_column(d, nms[i])) return false;
		} else if (x.itype[i] == 2) {
			std::vector<std::string> d = x.getS(i);
			if (!add_column(d, nms[i])) return false;
		} else if (x.itype[i] == 3) {
			std::vector<int8_t> d = x.getB(i);
			if (!add_column(d, nms[i])) return false;
		} else {
			std::vector<SpatTime_t> d = x.getT(i);
			if (!add_column(d, nms[i])) return false;
		}
	}
	return true;
}


bool SpatDataFrame::rbind(SpatDataFrame &x) {

	size_t nr1 = nrow();
	size_t nr2 = x.nrow();
//	size_t nc1 = ncol();
	size_t nc2 = x.ncol();

//first add new columns
	std::vector<std::string> nms = names;
	for (size_t i=0; i<nc2; i++) {
		int j = where_in_vector(x.names[i], nms, false);
		if (j < 0) { // not in df
			size_t b = x.iplace[i];
			add_column(x.itype[i], x.names[i]);
			if (x.itype[i] == 0) { 
				size_t a = dv.size()-1;
				dv[a].insert(dv[a].begin()+nr1, 
					x.dv[b].begin(), x.dv[b].end());
			} else if (x.itype[i] == 1) { 
				size_t a = iv.size()-1;
				iv[a].insert(iv[a].begin()+nr1, 
					x.iv[b].begin(), x.iv[b].end());
			} else if (x.itype[i] == 2) { 
				size_t a = sv.size()-1;
				sv[a].insert(sv[a].begin()+nr1, 
					x.sv[b].begin(), x.sv[b].end());
			} else if (x.itype[i] == 3) { 
				size_t a = bv.size()-1;
				bv[a].insert(bv[a].begin()+nr1, 
					x.bv[b].begin(), x.bv[b].end());
			} else { 
				size_t a = tv.size()-1;
				tv[a].insert(tv[a].begin()+nr1, 
					x.tv[b].begin(), x.tv[b].end());
			} 
		} else {
			size_t a = iplace[j];
			size_t b = x.iplace[i];
			if (itype[j] == x.itype[i]) {
				if (itype[j] == 0) { 
					dv[a].insert(dv[a].begin()+nr1, 
					x.dv[b].begin(), x.dv[b].end());
				} else if (itype[j] == 1) { 
					iv[a].insert(iv[a].begin()+nr1, 
					x.iv[b].begin(), x.iv[b].end());
				} else if (itype[j] == 2) { 
					sv[a].insert(sv[a].begin()+nr1, 
					x.sv[b].begin(), x.sv[b].end());
				} else if (itype[j] == 3) { 
					bv[a].insert(bv[a].begin()+nr1, 
					x.bv[b].begin(), x.bv[b].end());
				} else { 
					tv[a].insert(tv[a].begin()+nr1, 
					x.tv[b].begin(), x.tv[b].end());
				} 
			} else {
				if (itype[j] == 2) {
					sv[a].resize(nr1);
					if (x.itype[i] == 0) {
						for (size_t k=0; k<nr2; k++) {
							sv[a].push_back(std::to_string(x.dv[b][k]));
						}
					} else if (x.itype[i] == 1){
						for (size_t k=0; k<nr2; k++) {
							sv[a].push_back(std::to_string(x.iv[b][k]));
						}
					} else if (x.itype[i] == 3){
						for (size_t k=0; k<nr2; k++) {
							if (x.bv[b][k] == 0) {
								sv[a].push_back("TRUE");
							} else if (x.bv[b][k] == 1) {
								sv[a].push_back("FALSE");
							} else {
								sv[a].push_back("NAS");								
							}
						}
					} else {
						for (size_t k=0; k<nr2; k++) {
							sv[a].push_back(std::to_string(x.tv[b][k]));
						}
					}
				} else if (itype[j] == 0) { 
					dv[a].resize(nr1);
					if (x.itype[i] == 1) {
						for (size_t k=0; k<nr2; k++) {
							dv[a].push_back(x.iv[b][k]);
						}
					} else if (x.itype[i] == 2) {
						// could try to_double instead
						dv[a].resize(nr1 + nr2);
					} else if (x.itype[i] == 3) {
						for (size_t k=0; k<nr2; k++) {
							if (x.iv[b][k] < 2) {
								dv[a].push_back(x.bv[b][k]);
							} else {
								dv[a].push_back(NAN);								
							}
						}
					} else {
						for (size_t k=0; k<nr2; k++) {
							dv[a].push_back(x.tv[b][k]);
						}
					}

				} else if (itype[j] == 1) { 
					iv[a].resize(nr1);
					if (x.itype[i] == 0) {
						long longNA = NA<long>::value;
						for (size_t k=0; k<nr2; k++) {
							if (std::isnan(x.dv[b][k])) {
								iv[a].push_back(longNA);								
							} else {
								iv[a].push_back(x.dv[b][k]);
							}
						}
					} else if (x.itype[i] == 3) {
						long longNA = NA<long>::value;
						for (size_t k=0; k<nr2; k++) {
							if (x.bv[b][k] > 1) {
								iv[a].push_back(longNA);								
							} else {
								iv[a].push_back(x.bv[b][k]);
							}
						}
					} else {
						iv[a].resize(nr1 + nr2);						
					}
				} else if (itype[j] == 3) {					
					if (x.itype[i] == 1) {
						for (size_t k=0; k<nr2; k++) {
							if (x.dv[b][k] == 0 || x.dv[b][k] == 1) {
								bv[a].push_back(x.dv[b][k]);								
							} else {
								bv[a].push_back(2);
							}
						}
					} else {
						bv[a].resize(nr1 + nr2);						
					}
				} else if (itype[j] == 4) {
					tv[a].resize(nr1 + nr2);						
				}
			}
		}
	}

	resize_rows(nr1 + nr2);
	return true;
}


std::vector<std::string> SpatDataFrame::get_names() {
	return names;
}


void SpatDataFrame::set_names(std::vector<std::string> nms){
	if (ncol() == nms.size()) {
        make_valid_names(nms);
        make_unique_names(nms);
		names = nms;
	} else {
		setError("number of names is not correct");
	}
}


std::vector<std::string> SpatDataFrame::get_datatypes() {
	std::vector<std::string> types = {"double", "long", "string", "bool", "time"};
	std::vector<std::string> stype(itype.size());
	for (size_t i=0; i<itype.size(); i++) {
		stype[i] = types[itype[i]]; 
	}
	return stype;
}


std::string SpatDataFrame::get_datatype(std::string field) {
	int i = where_in_vector(field, get_names(), false);
	if (i < 0) return "";
	i = itype[i]; 
	std::vector<std::string> types = {"double", "long", "string", "bool", "time"};
	return types[i]; 
}

std::string SpatDataFrame::get_datatype(int field) {
	if ((field < 0) || (field > (int)(ncol()-1))) return "";
	std::vector<std::string> types = {"double", "long", "string", "bool", "time"};
	return types[itype[field]]; 
}

bool SpatDataFrame::field_exists(std::string field) {
	return is_in_vector(field, get_names()); 
}


int SpatDataFrame::get_fieldindex(std::string field) {
	return where_in_vector(field, get_names(), false);
}


// only doing this for one column for now
SpatDataFrame SpatDataFrame::unique(int col) {
	SpatDataFrame out = subset_cols(col);

	if (out.hasError()) return out;
	if (out.itype[0] == 0) {
		size_t sz = nrow();
		out.dv[0].erase(std::remove_if(out.dv[0].begin(), out.dv[0].end(),
            [](const double& value) { return std::isnan(value); }), out.dv[0].end());
		bool hasNAN = sz > out.dv[0].size();
		std::sort(out.dv[0].begin(), out.dv[0].end());
		out.dv[0].erase(std::unique(out.dv[0].begin(), out.dv[0].end()), out.dv[0].end());
		if (hasNAN) out.dv[0].push_back(NAN);

	} else if (out.itype[0] == 1) {
		std::sort(out.iv[0].begin(), out.iv[0].end());
		out.iv[0].erase(std::unique(out.iv[0].begin(), out.iv[0].end()), out.iv[0].end());
	} else if (out.itype[0] == 2) {
		std::sort(out.sv[0].begin(), out.sv[0].end());
		out.sv[0].erase(std::unique(out.sv[0].begin(), out.sv[0].end()), out.sv[0].end());
	} else if (out.itype[0] == 3) {
		std::sort(out.bv[0].begin(), out.bv[0].end());
		out.bv[0].erase(std::unique(out.bv[0].begin(), out.bv[0].end()), out.bv[0].end());
	} else if (out.itype[0] == 4) {
		std::sort(out.tv[0].begin(), out.tv[0].end());
		out.tv[0].erase(std::unique(out.tv[0].begin(), out.tv[0].end()), out.tv[0].end());
	}
	return out;
}


std::vector<int> SpatDataFrame::getIndex(int col, SpatDataFrame &x) {
	size_t nd = nrow();
	x = unique(col);
	size_t nu = x.nrow();
	std::vector<int> idx(nd, -1);
	size_t ccol = iplace[col];
	if (x.itype[0] == 0) {
		for (size_t i=0; i<nd; i++) {
			for (size_t j=0; j<nu; j++) {
				if (std::isnan(dv[ccol][i]) && (std::isnan(x.dv[0][j]))) {
					idx[i] = j;
					continue;
				} else if (dv[ccol][i] == x.dv[0][j]) {
					idx[i] = j;
					continue;
				}
			}
		}
	} else if (x.itype[0] == 1) {
		for (size_t i=0; i<nd; i++) {
			for (size_t j=0; j<nu; j++) {
				if (iv[ccol][i] == x.iv[0][j]) {
					idx[i] = j;
					continue;
				}
			}
		}
	} else if (x.itype[0] == 2) {
		for (size_t i=0; i<nd; i++) {
			for (size_t j=0; j<nu; j++) {
				if (sv[ccol][i] == x.sv[0][j]) {
					idx[i] = j;
					continue;
				}
			}
		}
	} else if (x.itype[0] == 3) {
		for (size_t i=0; i<nd; i++) {
			for (size_t j=0; j<nu; j++) {
				if (bv[ccol][i] == x.bv[0][j]) {
					idx[i] = j;
					continue;
				}
			}
		}
	} else if (x.itype[0] == 4) {
		for (size_t i=0; i<nd; i++) {
			for (size_t j=0; j<nu; j++) {
				if (tv[ccol][i] == x.tv[0][j]) {
					idx[i] = j;
					continue;
				}
			}
		}
	}
	return idx;
}


std::vector<double> SpatDataFrame::as_double(size_t v) {
	std::vector<double> out;
	if (v >= ncol()) {
		setError("attempting to read a column that does not exist");
		return out;
	}
	if (itype[v] == 2) {
		setError("as_double not available for string");
		return out;
	}
	size_t j = iplace[v];
	if (itype[v] == 0) return dv[j];
//	if (itype[v] == 1) {
	out.reserve(nrow());
	if (itype[v]==1) {
		long longNA = NA<long>::value;		
		for (size_t i=0; i<nrow(); i++){
			if (tv[j][i] == longNA) {
				out[i] = NAN;	
			} else {
				out[i] = (double)iv[j][i];			
			}
		}	
	} else if (itype[v]==3) {
		for (size_t i=0; i<nrow(); i++){
			if (bv[j][i] > 1) {
				out[i] = NAN;
			} else {
				out[i] = (double)bv[j][i];			
			}
		}
	} else if (itype[v]==4) {
		SpatTime_t timeNA = NA<SpatTime_t>::value;
		for (size_t i=0; i<nrow(); i++){
			if (tv[j][i] == timeNA) {
				out[i] = NAN;				
			} else {
				out[i] = (double)tv[j][i];
			}
		}
	} 
	return out;
}	
	

std::vector<long> SpatDataFrame::as_long(size_t v) {
	std::vector<long> out;
	if (v >= ncol()) {
		setError("attempting to read a column that does not exist");
		return out;
	}
	if (itype[v] == 2) {
		setError("as_long not available for string");
		return out;
	}
	size_t j = iplace[v];
	if (itype[v] == 1) return iv[j];
//	if (itype[v] == 0) {
	out.reserve(nrow());
	long longNA = NA<long>::value;		
	if (itype[v] == 0) {
		for (size_t i=0; i<nrow(); i++){
			if (std::isnan(dv[j][i])) {
				out.push_back(longNA);			
			} else {
				out.push_back( (long)dv[j][i] );
			}
		}
	} else if (itype[v]==3) {
		for (size_t i=0; i<nrow(); i++){
			if (bv[j][i] > 1) {
				out[i] = longNA;
			} else {
				out[i] = (long) bv[j][i];			
			}
		}
	} else if (itype[v]==4) {
		SpatTime_t timeNA = NA<SpatTime_t>::value;
		for (size_t i=0; i<nrow(); i++){
			if (tv[j][i] == timeNA) {
				out[i] = longNA;				
			} else {
				out[i] = (long) tv[j][i];
			}
		}
	} 
	return out;
}	
	

std::vector<std::string> SpatDataFrame::as_string(size_t v) {
	std::vector<std::string> out;
	if (v >= ncol()) {
		setError("attempting to read a column that does not exist");
		return out;
	}
	std::string dt = get_datatype(v);
	size_t j = iplace[v];

	if (dt == "string") return sv[j];
	out.reserve(nrow());
	if (dt == "double") {
		for (size_t i=0; i<nrow(); i++) {
			out.push_back(double_to_string(dv[j][i]));
		}
	} else if (dt == "long") {
		for (size_t i=0; i<nrow(); i++) {
			out.push_back(std::to_string(iv[j][i]));
		}
	}
	return out;
}

