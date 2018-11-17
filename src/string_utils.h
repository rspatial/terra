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

bool file_exists(const std::string& name);
std::string concatenate(std::vector<std::string> v, std::string delim);
void lowercase(std::string &s);
bool is_in_set(std::string s, std::vector<std::string> ss);
std::string is_in_set_default(std::string s, std::vector<std::string> ss, std::string defvalue, bool tolower);
std::vector<std::string> strsplit(std::string s, std::string delimiter);
std::vector<double> str2dbl(std::vector<std::string> s);
std::vector<std::string> dbl2str(std::vector<double> d);
std::string getFileExt(const std::string& s);
std::string setFileExt(const std::string& s, const std::string& ext);
std::string basename(std::string filename);
void lrtrim(std::string &s);
bool in_string(const std::string &x, std::string part);
