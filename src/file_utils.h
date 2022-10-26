// Copyright (c) 2018-2022  Robert J. Hijmans
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

#ifndef FILEUTILS_GUARD
#define FILEUTILS_GUARD

bool file_exists(const std::string& name);
bool filepath_exists(const std::string& name);
bool can_write(std::vector<std::string> filenames, std::vector<std::string> srcnames, bool overwrite, std::string &msg);
std::string getFileExt(const std::string& s);
std::string setFileExt(const std::string& s, const std::string& ext);
std::string basename(std::string filename);
std::string basename_noext(std::string filename);
std::string noext(std::string filename);
std::string tempFile(std::string tmpdir, unsigned pid, std::string ext);
std::string dirname(std::string filename);
bool write_text(std::string filename, std::vector<std::string> s);
std::vector<std::string> read_text(std::string filename);

#endif

