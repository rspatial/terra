// Copyright (c) 2018-2020  Robert J. Hijmans
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

#ifndef SPATMESSAGESGUARD
#define SPATMESSAGESGUARD

#include <string>

class SpatMessages {
	public:
		bool success = true;
		bool has_error = false;
		bool has_warning = false;
		std::string error;
		std::vector<std::string> warnings;

		void setError(std::string s) {
			has_error = true;
			success = false;
			error = s;
		}
		std::string getError() {
			has_error = false;
			success = true;
			return error;
		}
		
		void addWarning(std::string s) {
			has_warning = true;
			warnings.push_back(s);
		}
		
		std::vector<std::string> getMessages() {
			std::string warns = "";
			if (warnings.size() > 0) {
				warns = warnings[0];
				for (size_t i = 1; i<warnings.size(); i++) {
					warns = warns  + "\n" + warnings[i];
				}
			}
			std::vector<std::string> msg = { error, warns};
			has_error = false;
			has_warning = false;
			error = "";
			warnings.resize(0);
			return msg;
		}
};

#endif