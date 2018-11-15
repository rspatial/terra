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

#include "spatBase.h"

SpatOptions::SpatOptions() {}


SpatOptions::SpatOptions(const SpatOptions &opt) {
	data_type = opt.data_type;
	file_type = opt.data_type; 
	file_name = "";
	temp_dir = opt.temp_dir;
	mem_frac = opt.mem_frac;
}


SpatOptions SpatOptions::deepcopy(const SpatOptions &opt) {
	return SpatOptions(opt);
}

std::string SpatOptions::get_datatype() {
	return data_type;
}

std::string SpatOptions::get_filetype() {
	return file_type;
}

bool SpatOptions::get_overwrite() { return over_write; }

		
// check if valid datatype`	
void SpatOptions::set_overwrite(bool d) { over_write = d; }

void SpatOptions::set_datatype(std::string d) { data_type = d; }

void SpatOptions::set_filetype(std::string d) { file_type = d; }

void SpatOptions::set_filename(std::string d) { file_name = d; }

void SpatOptions::set_tempdir(std::string d) {
	// check if exists
	temp_dir = d;
}


std::string SpatOptions::get_tempdir() { return temp_dir; }
std::string SpatOptions::get_filename() { return file_name; }

void SpatOptions::set_memfrac(double d) {
	if ((d >= 0.1) && (d < 0.9)) { 
		mem_frac = d;
		return;
	} 
	//setError;
}

double SpatOptions::get_memfrac() { return mem_frac; }
		