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

