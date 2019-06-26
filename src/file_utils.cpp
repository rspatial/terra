#include "SpatBase.h"
#include <fstream>



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




bool file_exists(const std::string& name) {
	std::ifstream f(name.c_str());
	return f.good();
}


bool canWrite(std::string filename) {
	FILE *fp = fopen(filename.c_str(), "w");
	if (fp == NULL) {
		return false;
	}
	fclose(fp);
	return true;
}

SpatMessages can_write(std::string filename, bool overwrite) {
	SpatMessages msg;
	if (file_exists(filename)) {
		if (overwrite) {
			remove(filename.c_str());
		} else {
			msg.setError("file exists");
			return msg;
		}
	}
	if (!canWrite(filename)) {
		msg.setError("cannot write file");
	}
	return msg;
}


