#include <vector>
#include <fstream>


bool writeINT2(std::string filename, std::vector<double> v) {
	std::vector<short> values(v.begin(), v.end());
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(short));
	file.close();
	return(true);
}

bool writeINT4(std::string filename, std::vector<double> v) {
	std::vector<long> values(v.begin(), v.end());
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(long));
	file.close();
	return(true);
}

bool writeFLT4(std::string filename, std::vector<double> v) {
	std::vector<float> values(v.begin(), v.end());
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(float));
	file.close();
	return(true);
}

bool writeFLT8(std::string filename, std::vector<double> v) {
	std::vector<double> values(v.begin(), v.end());
	std::ofstream file(filename, std::ios::out | std::ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(double));
	file.close();
	return(true);
}

