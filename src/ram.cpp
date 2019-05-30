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

#ifdef _WIN32
#include <windows.h>
#elif __linux__
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <unistd.h>
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
		struct sysinfo memInfo;
		sysinfo (&memInfo);
		ram = memInfo.freeram;
	#elif __APPLE__
		vm_size_t page_size;
		mach_port_t mach_port;
		mach_msg_type_number_t count;
		vm_statistics64_data_t vm_stats;
		mach_port = mach_host_self();
		count = sizeof(vm_stats) / sizeof(natural_t);
		if (KERN_SUCCESS == host_page_size(mach_port, &page_size) &&
					KERN_SUCCESS == host_statistics64(mach_port, HOST_VM_INFO,
											(host_info64_t)&vm_stats, &count)) {
			ram = (int64_t)vm_stats.free_count * (int64_t)page_size;
		}
	#endif
	return ram / 8;
}


