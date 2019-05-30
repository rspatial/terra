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

#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <fstream>



bool file_exists(const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}

std::string concatenate(std::vector<std::string> v, std::string delim) {
	for (size_t i=0; i<(v.size()-1); i++) {
		v[i] = v[i] + delim;
	}
	std::string s;
	for (const auto &piece : v) s += piece;
	return s;
}

bool in_string(const std::string &x, std::string part) {
	size_t f = x.find(part);
	if (f == std::string::npos) {
		return false;
	} else {
		return true;
	}
}


void lowercase(std::string &s) {
	std::transform(s.begin(), s.end(), s.begin(), ::tolower);
}

bool is_in_set(std::string s, std::vector<std::string> ss) {
	std::set<std::string> sset (ss.begin(), ss.end());
	return sset.find(s) != sset.end();
}

std::string is_in_set_default(std::string s, std::vector<std::string> ss, std::string defvalue, bool tolower) {
	if (tolower) lowercase(s);
	std::set<std::string> sset (ss.begin(), ss.end());
	if (sset.find(s) == sset.end() ) {
		s = defvalue;
	}
	return s;
}


std::vector<std::string> strsplit(std::string s, std::string delimiter){
	std::vector<std::string> out;
	size_t pos = 0;
	std::string token;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		out.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	token = s.substr(0, pos);
	out.push_back(token);
	return out;
}


std::vector<double> str2dbl(std::vector<std::string> s) {
	std::vector<double> d (s.size());
	std::transform(s.begin(), s.end(), d.begin(), [](const std::string& val) {
		return std::stod(val);
	});
	return d;
}

std::vector<std::string> dbl2str(std::vector<double> d) {
	std::vector<std::string> s (d.size());
	std::transform(d.begin(), d.end(), s.begin(),
			[](double i) { return std::to_string(i); }
	);
	return s;
}



std::string getFileExt(const std::string& s) {
	size_t i = s.rfind('.', s.length());
	if (i != std::string::npos) {
		return(s.substr(i, s.length() - i));
	}
	return("");
}

std::string setFileExt(const std::string& s, const std::string& ext) {
	size_t i = s.rfind('.', s.length());
	if (i != std::string::npos) {
		return(s.substr(0, i) + ext);
	}
	return(s + ext);
}

std::string basename(std::string filename) {
	const size_t i = filename.find_last_of("\\/");
	if (std::string::npos != i) {
		filename.erase(0, i + 1);
	}
	const size_t p = filename.rfind('.');
	if (std::string::npos != p) {
		filename.erase(p);
	}
	return filename;
}



// trim from start (in place)
void ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

// trim from both ends (in place)
void lrtrim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

// trim from start (copying)
std::string ltrim_copy(std::string s) {
    ltrim(s);
    return s;
}

// trim from end (copying)
std::string rtrim_copy(std::string s) {
    rtrim(s);
    return s;
}

// trim from both ends (copying)
std::string lrtrim_copy(std::string s) {
    lrtrim(s);
    return s;
}



