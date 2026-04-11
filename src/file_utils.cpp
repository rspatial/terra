// Copyright (c) 2018-2026  Robert J. Hijmans
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
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>
#include <thread>

#include <sys/types.h>
#include <sys/stat.h>

#include "cpl_vsi.h"

static inline bool is_vsi(const std::string& path) {
	return path.size() > 4 && path.substr(0, 4) == "/vsi";
}

/*
#if defined __has_include
#	if __has_include (<filesystem>)
# 		include <filesystem>
		namespace filesyst = std::filesystem;
#	else
#		include <experimental/filesystem>
		namespace filesyst = std::experimental::filesystem;
#	endif
#elif defined __GNUC__
#	if __GNUC__ < 8
#		include <experimental/filesystem>
		namespace filesyst = std::experimental::filesystem;
#	else 
# 		include <filesystem>
		namespace filesyst = std::filesystem;	
#	endif
#else 
#	include <filesystem>
    namespace filesyst = std::filesystem;
#endif
*/

bool write_text(std::string filename, std::vector<std::string> s) {
	std::ofstream f;
	f.open(filename);
	if (f.is_open()) {
		for (size_t i=0; i<s.size(); i++) {
			f << s[i] << std::endl;
		}
		f.close();
		return true;
	} else {
		return false;
	}
}


std::vector<std::string> read_text(std::string filename) {
	std::vector<std::string> s;
	if (is_vsi(filename)) {
		VSILFILE *fp = VSIFOpenL(filename.c_str(), "r");
		if (fp != nullptr) {
			char buf[4096];
			std::string residual;
			while (true) {
				size_t nread = VSIFReadL(buf, 1, sizeof(buf), fp);
				if (nread == 0) break;
				residual.append(buf, nread);
			}
			VSIFCloseL(fp);
			std::istringstream iss(residual);
			std::string line;
			while (std::getline(iss, line)) {
				s.push_back(line);
			}
		}
		return s;
	}
	std::string line;
	std::ifstream f(filename);
	if (f.is_open())  {
		while (getline(f, line)) {
			if (line.empty()) {
				s.push_back("");
			} else {
				s.push_back(line);
			}
		}
		f.close();
	}
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

std::string noext(std::string filename) {
	const size_t p = filename.rfind('.');
	if (std::string::npos != p) {
		filename.erase(p);
	}
	return filename;
}

std::string basename(std::string filename) {
	const size_t i = filename.find_last_of("\\/");
	if (std::string::npos != i) {
		filename.erase(0, i + 1);
	}
	return filename;
}


std::string basename_noext(std::string filename) {
	filename = basename(filename);
	filename = noext(filename);
	return filename;
}


std::string dirname(std::string filename) {
	const size_t i = filename.find_last_of("\\/");
	if (std::string::npos != i) {
		return( filename.substr(0, i) );
	} else {
		return ("");
	}
}

// For /vsizip/, /vsigzip, extract the path in the .zip or .gz file.
std::string get_vsi_container(const std::string& path) {
	std::string prefix;
	if (path.size() > 8 && path.substr(0, 8) == "/vsizip/") {
		prefix = "/vsizip/";
	} else if (path.size() > 10 && path.substr(0, 10) == "/vsigzip/") {
		prefix = "/vsigzip/";
	} else {
		return "";
	}
	std::string rest = path.substr(prefix.size());
	// Curly brace format: /vsizip/{container_path}/inner
	if (!rest.empty() && rest[0] == '{') {
		size_t end = rest.find('}');
		if (end != std::string::npos) {
			return rest.substr(1, end - 1);
		}
	}
	// Standard format: look for .zip or .gz extension
	for (const auto& ext : {".zip", ".ZIP", ".gz", ".GZ"}) {
		size_t pos = rest.find(ext);
		if (pos != std::string::npos) {
			return rest.substr(0, pos + strlen(ext));
		}
	}
	return "";
}


bool file_exists(const std::string& name) {
	if (is_vsi(name)) {
		std::string container = get_vsi_container(name);
		if (!container.empty()) {
			std::ifstream cf(container.c_str());
			if (!cf.good()) return false;
		}
		VSIStatBufL statBuf;
		return VSIStatL(name.c_str(), &statBuf) == 0;
	}
	std::ifstream f(name.c_str());
	return f.good();
}


bool path_exists(std::string path) {

/*
	filesyst::path filepath = path;
	return filesyst::exists(filepath);
*/
	struct stat info;
	stat(path.c_str(), &info);
	if (info.st_mode & S_IFDIR) {
		return true;
	}
	return false;
}



bool canWrite(std::string filename) {
	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == NULL) {
		return false;
	}
	fclose(fp);
	remove(filename.c_str());
	return true;
}


std::string get_path(const std::string filename) {
	size_t found = filename.find_last_of("/\\");
	std::string result = filename.substr(0, found);
	return result;
}

bool filepath_exists(const std::string& name) {
	std::string p = get_path(name);
	return path_exists(p);
}




/*
# c++17
#include <experimental/filesystem>
bool SpatRaster::differentFilenames(std::vector<std::string> outf) {
	std::vector<std::string> inf = filenames();
	for (size_t i=0; i<inf.size(); i++) {
		if (inf[i] == "") continue;
		std::experimental::filesystem::path pin = inf[i];
		for (size_t j=0; j<outf.size(); j++) {
			std::experimental::filesystem::path pout = outf[i];
			if (pin.compare(pout) == 0) return false;
		}
	}
	return true;
}
*/

bool differentFilenames(std::vector<std::string> inf, std::vector<std::string> outf, std::string &msg) {
	#ifdef _WIN32
	for (size_t j=0; j<outf.size(); j++) {
		std::transform(outf[j].begin(), outf[j].end(), outf[j].begin(), ::tolower);
	}
	#endif

	for (size_t i=0; i<inf.size(); i++) {
		if (inf[i].empty()) continue;
		#ifdef _WIN32
		std::transform(inf[i].begin(), inf[i].end(), inf[i].begin(), ::tolower);
		#endif
		for (size_t j=0; j<outf.size(); j++) {
			if (inf[i] == outf[j]) {
				msg = "source and target filename cannot be the same";
				return false;
			}
		}
	}
	size_t n = outf.size();
	std::sort( outf.begin(), outf.end() );
	outf.erase(std::unique(outf.begin(), outf.end()), outf.end());
	if (n > outf.size()) {
		msg = "duplicate filenames";
		return false;
	}
	return true;
}


bool can_write(std::vector<std::string> filenames, std::vector<std::string> srcnames, bool overwrite, std::string &msg) {

	if (!differentFilenames(srcnames, filenames, msg)) {
		return false;
	}

	for (size_t i=0; i<filenames.size(); i++) {
		if (!filenames[i].empty() && file_exists(filenames[i])) {
			if (overwrite) {
				bool removed = false;
				if (is_vsi(filenames[i])) {
					// for archive VSI paths, delete the container file
					std::string container = get_vsi_container(filenames[i]);
					if (!container.empty()) {
						removed = (remove(container.c_str()) == 0);
					}
					if (!removed) {
						removed = (VSIUnlink(filenames[i].c_str()) == 0);
					}
				} else {
					removed = (remove(filenames[i].c_str()) == 0);
				}
				if (!removed) {
					msg = "cannot overwrite existing file";
					return false;
				}
				if (!is_vsi(filenames[i])) {
					std::vector<std::string> exts = {".vat.dbf", ".vat.cpg", ".json", ".aux.xml"};
					for (size_t j=0; j<exts.size(); j++) {
						std::string f = filenames[i] + exts[j];
						if (file_exists(f)) {
							remove(f.c_str());
						}
					}
				}
			} else {
				msg = "file exists. You can use 'overwrite=TRUE' to overwrite it";
				return false;
			}
		} else if (!canWrite(filenames[i])) {
			if (filenames[i].substr(0, 4) == "/vsi") continue; 
			std::string path = get_path(filenames[i]);
			if (!path_exists(path)) {
				msg = "path does not exist";
			} else {
				msg = "cannot write file";
			}
			return false;
		}
	}
	return true;
}

#include "common.h"

std::mt19937 my_rgen;

std::string tempFile(std::string tmpdir, std::string fname, std::string ext) {


    std::vector<char> characters = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K',
    'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m',
    'n','o','p','q','r','s','t','u','v','w','x','y','z' };

    std::uniform_int_distribution<std::mt19937::result_type> rand_nr(0, characters.size()-1); 
	
    std::string randname;
	randname.reserve(15);
    for (int i = 0; i < 15; i++) {
        randname += characters[rand_nr(my_rgen)];
    }
  
	std::string filename =  tmpdir + "/spat_" + fname + "_" + randname + ext;
	if (file_exists(filename)) {
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
		filename = tempFile(tmpdir, fname, ext);
	}
	return filename;
}



/*
std::string tempFile(std::string tmpdir, unsigned pid, std::string ext) {
    std::vector<char> characters = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K',
    'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m',
    'n','o','p','q','r','s','t','u','v','w','x','y','z' };
    std::uniform_int_distribution<> distrib(0, characters.size()-1);
    auto draw = [ characters, &distrib, &generator ]() {
		return characters[ distrib(generator) ];
	};
    std::string filename(15, 0);
    std::generate_n(filename.begin(), 15, draw);
	filename = tmpdir + "/spat_" + filename + "_" + std::to_string(pid) + ext;
	if (file_exists(filename)) {
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
		return tempFile(tmpdir, pid, ext);
	}
	return filename;
}
*/


// --- open-file limit --------------------------------------------------------
// Returns three numbers describing the file-descriptor situation:
//   nopen  – how many handles/descriptors the process currently has open
//   soft   – the current (soft) per-process limit
//   hard   – the maximum the soft limit can be raised to (without privileges)
//
// POSIX : limits from getrlimit(RLIMIT_NOFILE), count from /proc/self/fd
//         (Linux) or /dev/fd (macOS).
// Windows: GDAL uses Win32 CreateFile (kernel handles), not CRT fopen, so we
//          count kernel handles with GetProcessHandleCount.  The per-process
//          handle limit on Windows is very large (16M+), so in practice the
//          readStart guard will only trigger on POSIX systems.

#ifdef _WIN32
// windows.h already included at the top of ram.cpp; here we may need it too
#ifndef _WINDOWS_
#include <windows.h>
#endif
#else
#include <sys/resource.h> // getrlimit, RLIMIT_NOFILE
#include <dirent.h>       // opendir, readdir, closedir
#endif

void open_file_limit(size_t &nopen, size_t &soft, size_t &hard) {

#ifdef _WIN32
	DWORD hcount = 0;
	GetProcessHandleCount(GetCurrentProcess(), &hcount);
	nopen = (size_t) hcount;
	soft  = 16777216;
	hard  = 16777216;
#else
	struct rlimit rl;
	if (getrlimit(RLIMIT_NOFILE, &rl) == 0) {
		soft = (size_t) rl.rlim_cur;
		hard = (rl.rlim_max == RLIM_INFINITY) ? (size_t) -1 : (size_t) rl.rlim_max;
	} else {
		soft = 256;
		hard = 256;
	}

	nopen = 0;
#ifdef __linux__
	const char *fddir = "/proc/self/fd";
#else
	const char *fddir = "/dev/fd";
#endif
	DIR *d = opendir(fddir);
	if (d) {
		while (readdir(d)) nopen++;
		closedir(d);
		if (nopen >= 3) nopen -= 3;  // subtract ".", "..", and the opendir() fd
	}
#endif
}

