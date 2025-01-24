// Copyright (c) 2018-2025  Robert J. Hijmans
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

#ifdef _WIN32
#include <windows.h>
#elif __linux__
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#elif __APPLE__
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#endif


double availableRAM() {
//https://stackoverflow.com/questions/38490320/how-to-query-amount-of-allocated-memory-on-linux-and-osx
//https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
	double ram = 1e+9;
	// return available RAM in number of double (8 byte) cells.
	#ifdef _WIN32
		MEMORYSTATUSEX statex;
		statex.dwLength = sizeof(statex);
		GlobalMemoryStatusEx(&statex);
		ram = statex.ullAvailPhys;
	#elif __linux__
		unsigned long memAvailable = 0;
		std::ifstream meminfo("/proc/meminfo");
		std::string line;
		while (std::getline(meminfo, line)) {
			std::istringstream iss(line);
			std::string key;
			if (std::getline(iss, key, ':')) {
				if (key == "MemAvailable") {
					std::string value_str;
					if (std::getline(iss >> std::ws, value_str)) {
						std::stringstream value_stream(value_str);
						unsigned long mem_available_kb;
						std::string units;
						if (value_stream >> mem_available_kb >> units) {
							memAvailable = mem_available_kb;
						}
					}
					break;
				}
			}

		}
		if (memAvailable == 0) {
			struct sysinfo memInfo;
			sysinfo (&memInfo);
			ram = memInfo.freeram;
		} else {
			ram = memAvailable * 1024;
		}
	#elif __APPLE__

		vm_size_t page_size;
		mach_port_t mach_port;
		mach_msg_type_number_t count;
	#if defined(__ppc__) || defined(__i386__)
		vm_statistics_data_t vm_stats;
	#else
		vm_statistics64_data_t vm_stats;
	#endif

		mach_port = mach_host_self();
		count = sizeof(vm_stats) / sizeof(natural_t);
	#if defined(__ppc__) || defined(__i386__)
		if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
			KERN_SUCCESS == host_statistics(mach_port, HOST_VM_INFO,
											(host_info_t)&vm_stats, &count)) {
			long long free_memory = ((int32_t)vm_stats.free_count +
                               (int32_t)vm_stats.inactive_count) * (int32_t)page_size;
	#else
		if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
			KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
											(host_info64_t)&vm_stats, &count)) {
			long long free_memory = ((int64_t)vm_stats.free_count +
                               (int64_t)vm_stats.inactive_count) * (int64_t)page_size;
	#endif
			ram = free_memory;
		//https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
		}

	#endif
	return ram / 8;  // 8 bytes for each double
}
