using namespace std;

#include <vector>
#include <fstream>


bool writeINT2(string filename, std::vector<double> v) {
	std::vector<short> values(v.begin(), v.end());
	std::ofstream file(filename, ios::out | ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(short));
	file.close();
	return(true);
}

bool writeINT4(string filename, std::vector<double> v) {
	std::vector<long> values(v.begin(), v.end());
	std::ofstream file(filename, ios::out | ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(long));
	file.close();
	return(true);
}

bool writeFLT4(string filename, std::vector<double> v) {
	std::vector<float> values(v.begin(), v.end());
	std::ofstream file(filename, ios::out | ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(float));
	file.close();
	return(true);
}

bool writeFLT8(string filename, std::vector<double> v) {
	std::vector<double> values(v.begin(), v.end());
	std::ofstream file(filename, ios::out | ios::binary);
	file.write((char*)&values[0], values.size() * sizeof(double));
	file.close();
	return(true);
}

