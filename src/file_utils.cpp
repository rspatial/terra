#include "SpatBase.h"
#include <fstream>
#include <random>
#include <chrono>



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


std::string tempFile(std::string tmpdir, std::string ext) {
    std::vector<char> characters = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K',
    'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m',
    'n','o','p','q','r','s','t','u','v','w','x','y','z' };
    std::default_random_engine generator(std::random_device{}());
	double seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed(seed);
    std::uniform_int_distribution<> distrib(0, characters.size()-1);
    auto draw = [ characters, &distrib, &generator ]() {
		return characters[ distrib(generator) ];
	};
    std::string filename(15, 0);
    std::generate_n(filename.begin(), 15, draw);
	filename = tmpdir + "/spat_" + filename + ext;
	return filename;
}

