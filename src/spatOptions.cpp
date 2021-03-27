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

#include "spatBase.h"
#include "string_utils.h"


SpatOptions::SpatOptions() {}


SpatOptions::SpatOptions(const SpatOptions &opt) {
	tempdir = opt.tempdir;
	memfrac = opt.memfrac;
	todisk = opt.todisk;
	def_datatype = opt.def_datatype;
	def_filetype = opt.def_filetype; 
	filenames = {""};
	overwrite = false;
	progress = opt.progress;
	ncopies = opt.ncopies;
	verbose = opt.verbose;
	statistics = opt.statistics;
	steps = opt.steps;
	minrows = opt.minrows;
	//ncdfcopy = opt.ncdfcopy;
}

SpatOptions SpatOptions::deepCopy() {
	SpatOptions opt = *this;
	return opt;
}


//SpatOptions SpatOptions::deepCopy(const SpatOptions &opt) {
//	return SpatOptions(opt);
//}

//void SpatOptions::set_def_bandorder(std::string d) { def_bandorder = d; }
//std::string SpatOptions::get_def_bandorder() { return def_bandorder; }
//void SpatOptions::set_bandorder(std::string d) { bandorder = d; }
//std::string SpatOptions::get_bandorder() {if (bandorder != "") {return bandorder;} else {return def_datatype;}}

void SpatOptions::set_def_datatype(std::string d) { 
	std::vector<std::string> ss = {"INT1U", "INT2U", "INT4U", "INT2S", "INT4S", "FLT4S", "FLT8S" } ;
	if (is_in_vector(d, ss)) def_datatype = d; 
}
std::string SpatOptions::get_def_datatype() { return def_datatype; }

void SpatOptions::set_datatype(std::string d) { 
	std::vector<std::string> ss = {"INT1U", "INT2U", "INT4U", "INT2S", "INT4S", "FLT4S", "FLT8S" };
	if (is_in_vector(d, ss)) datatype = d; 
}
std::string SpatOptions::get_datatype() {if (datatype != "") {return datatype;} else {return def_datatype;}}

void SpatOptions::set_def_filetype(std::string d) { def_filetype = d; }
std::string SpatOptions::get_def_filetype() { return def_filetype;}

void SpatOptions::set_filetype(std::string d) { filetype = d; }
std::string SpatOptions::get_filetype() { return filetype;}

bool SpatOptions::get_overwrite() { return overwrite; }
void SpatOptions::set_overwrite(bool b) { overwrite = b; }

int SpatOptions::get_statistics() { return statistics; }
void SpatOptions::set_statistics(int s) { if ((s> 0) && (s<7)) statistics = s; }

//bool SpatOptions::get_ncdfcopy() { return ncdfcopy;}
//void SpatOptions::set_ncdfcopy(bool x) { ncdfcopy = x; }

void SpatOptions::set_def_verbose(bool v) { def_verbose = v; }
bool SpatOptions::get_def_verbose() { return def_verbose; }
bool SpatOptions::get_verbose() { return verbose; }
void SpatOptions::set_verbose(bool v) { verbose = v; }

bool SpatOptions::has_NAflag(double &flag) { 
	flag = NAflag;
	return hasNAflag; 
}

double SpatOptions::get_NAflag() { 
	return NAflag;
}

void SpatOptions::set_NAflag(double flag) { 
	NAflag = flag; 
	hasNAflag = true;
}

unsigned SpatOptions::get_progress() { return progress; }
void SpatOptions::set_progress(unsigned p) { 
	progress = p; 
}

bool SpatOptions::show_progress(unsigned n) { 
	return ((progress > 0) & (progress <= n));
}


//void SpatOptions::set_filename(std::string f) { 
//	f = lrtrim_copy(f); 
//	filenames = {f}; 
//}

void SpatOptions::set_filenames(std::vector<std::string> f) { 
	for (size_t i=0; i<f.size(); i++) {
		f[i] = lrtrim_copy(f[i]); 
	}
	filenames = f; 
}

std::string SpatOptions::get_filename() { 
	if (!filenames.empty() ) {
		return filenames[0]; 
	} else {
		return "";
	}
}


std::vector<std::string> SpatOptions::get_filenames() { 
	if (!filenames.empty() ) {
		return filenames; 
	} else {
		return {""};
	}
}

std::string SpatOptions::get_tempdir() { return tempdir; }

void SpatOptions::set_tempdir(std::string d) {
	// check if exists?
	tempdir = d;
}

double SpatOptions::get_memfrac() { return memfrac; }

void SpatOptions::set_memfrac(double d) {
	// allowing very high values for testing purposes
	if ((d >= 0.1) && (d <= 100)) { 
		memfrac = d;
		return;
	} 
	//setError;
}

bool SpatOptions::get_todisk() { return todisk; }
void SpatOptions::set_todisk(bool b) { todisk = b; }


void SpatOptions::set_steps(size_t n) { steps = std::max((size_t)1, n); }
size_t SpatOptions::get_steps(){ return steps; }

void SpatOptions::set_ncopies(size_t n) { ncopies = std::max((size_t)1, n); }
size_t SpatOptions::get_ncopies(){ return ncopies; }
