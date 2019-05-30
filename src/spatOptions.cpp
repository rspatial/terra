// Copyright (c) 2018-2019  Robert J. Hijmans
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
	filename = "";
	overwrite = false;	
	progress = opt.progress;
	blocksizemp = opt.blocksizemp;
}


SpatOptions SpatOptions::deepcopy(const SpatOptions &opt) {
	return SpatOptions(opt);
}

void SpatOptions::set_def_datatype(std::string d) { def_filetype = d; }
std::string SpatOptions::get_def_datatype() { return def_datatype; }

void SpatOptions::set_datatype(std::string d) { datatype = d; }
std::string SpatOptions::get_datatype() {if (datatype != "") {return datatype;} else {return def_datatype;}}

void SpatOptions::set_def_filetype(std::string d) { def_filetype = d; }
std::string SpatOptions::get_def_filetype() { return def_filetype;}
void SpatOptions::set_filetype(std::string d) { filetype = d; }
std::string SpatOptions::get_filetype() { return filetype;}

bool SpatOptions::get_overwrite() { return overwrite; }
void SpatOptions::set_overwrite(bool b) { overwrite = b; }

unsigned SpatOptions::get_progress() { return progress; }
void SpatOptions::set_progress(unsigned p) { 
	progress = p; 
}
bool SpatOptions::do_progress(unsigned n) { 
	return ((progress > 0) & (progress <= n));
}

unsigned SpatOptions::get_blocksizemp() { return blocksizemp; }
void SpatOptions::set_blocksizemp(unsigned x) { 
	blocksizemp = x; 
}


void SpatOptions::set_filename(std::string d) { lrtrim(d); filename = d; }
std::string SpatOptions::get_filename() { return filename; }

std::string SpatOptions::get_tempdir() { return tempdir; }
void SpatOptions::set_tempdir(std::string d) {
	// check if exists
	tempdir = d;
}

double SpatOptions::get_memfrac() { return memfrac; }
void SpatOptions::set_memfrac(double d) {
	if ((d >= 0.1) && (d < 0.9)) { 
		memfrac = d;
		return;
	} 
	//setError;
}

bool SpatOptions::get_todisk() { return todisk; }
void SpatOptions::set_todisk(bool b) { todisk = b; }
