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

#ifndef STRINGUTILS_GUARD
#define STRINGUTILS_GUARD

#include<string>
#include<vector>

void unquote(std::string &s);

std::string double_to_string(double x);
std::vector<std::string> double_to_string(const std::vector<double> &x, std::string prep);
std::vector<char *> string_to_charpnt(std::vector<std::string> s);

std::string concatenate(std::vector<std::string> v, std::string delim);
void lowercase(std::string &s);
void lowercase(std::vector<std::string> &ss);
std::string lower_case(std::string s);


std::string is_in_set_default(std::string s, std::vector<std::string> ss, std::string defvalue, bool tolower);
int where_in_set(std::string s, std::vector<std::string> ss, bool tolower);
bool is_in_vector(std::string s, std::vector<std::string> ss);
int where_in_vector(std::string s, const std::vector<std::string> &ss, const bool &tolower);

std::vector<std::string> getlastpart (std::vector<std::string> s, std::string delim);

std::vector<std::string> strsplit(std::string s, std::string delimiter);
std::vector<std::string> strsplit_first(std::string s, std::string delimiter);

std::vector<double> str2dbl(std::vector<std::string> s);
std::vector<int> str2int(std::vector<std::string> s);
std::vector<long> str2long(std::vector<std::string> s);

std::vector<std::string> dbl2str(std::vector<double> d);
void lrtrim(std::string &s);
std::string lrtrim_copy(std::string s);

bool in_string(const std::string &x, std::string part);
bool ends_on(std::string const &s, std::string const &end);

void make_unique_names(std::vector<std::string> &s);
void make_valid_names(std::vector<std::string> &s);

void str_replace(std::string& str, const std::string& from, const std::string& to);
size_t str_replace_all(std::string& str, const std::string& from, const std::string& to);


#endif

