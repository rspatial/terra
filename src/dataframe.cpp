#include "dataframe.h"
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include "util.h"


SpatDataFrame SpatDataFrame::subsetrows(std::vector<unsigned> range) { 
	SpatDataFrame out;
	out.names = names;
	out.itype = itype;
	out.iplace = iplace;
	out.dv.resize(dv.size());
	out.iv.resize(iv.size());
	out.sv.resize(sv.size());
	for (size_t i=0; i < range.size(); i++) {
		for (size_t j=0; j < dv.size(); j++) {
			out.dv[j].push_back(dv[j][i]);
		}
		for (size_t j=0; j < iv.size(); j++) {
			out.iv[j].push_back(iv[j][i]);
		}
		for (size_t j=0; j < sv.size(); j++) {
			out.sv[j].push_back(sv[j][i]);
		}
	}
	return out;	
}	

SpatDataFrame SpatDataFrame::subsetcols(std::vector<unsigned> range) { 
	SpatDataFrame out;
	out.dv.resize(0);
	out.iv.resize(0);
	out.sv.resize(0);
	unsigned dcnt=0;
	unsigned icnt=0;
	unsigned scnt=0;
	unsigned t, p;
	for (size_t i=0; i < range.size(); i++) {
		unsigned j = range[i];
		t = itype[j];
		p = iplace[j];
		out.names.push_back(names[j]);
		if (itype[j] == 0) {
			out.dv.push_back(dv[p]);
			out.iplace.push_back(dcnt);
			dcnt++;
		} else if (itype[j] == 1) {
			out.iv.push_back(iv[p]);
			out.iplace.push_back(icnt);
			icnt++;				
		} else {
			out.sv.push_back(sv[p]);
			out.iplace.push_back(scnt);
			scnt++;				
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
		} else {
			n = sv[0].size();			
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
}

	
	