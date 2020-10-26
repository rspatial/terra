// Robert Hijmans with improvements by Ben Fasoli
// https://github.com/rspatial/raster/pull/175

#ifdef _WIN32 
#include <windows.h>
#elif __linux__
#include <stdio.h>
#elif __APPLE__
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#endif

double availableRAM() {
	double ram = 1e+9;
	// return available RAM
	#ifdef _WIN32
		MEMORYSTATUSEX statex;
		statex.dwLength = sizeof(statex);
		GlobalMemoryStatusEx(&statex);
		ram = statex.ullAvailPhys;
	#elif __linux__
		// source available memory from /proc/meminfo
		// default to searching for MemAvailable field (kernel versions >= 3.14)
		FILE *fp = popen("awk '/MemAvailable/ {print $2}' /proc/meminfo", "r");
		if (fp == NULL) {
			return ram;
		}
		double ramkb;
		int ok = fscanf(fp, "%lf", &ramkb);  // returned in kB
		pclose(fp);
		if (ramkb > 0) {
			return ramkb * 1000.;
		}
		
		// fallback to estimating memory from other fields if MemAvailable not found
		FILE *fp2 = popen("awk -v low=$(grep low /proc/zoneinfo | awk '{k+=$2}END{print k}') '{a[$1]=$2}END{print a[\"MemFree:\"]+a[\"Active(file):\"]+a[\"Inactive(file):\"]+a[\"SReclaimable:\"]-(12*low);}' /proc/meminfo", "r");
		if (fp2 == NULL) {
			return ram;
		}
		ok = fscanf(fp2, "%lf", &ramkb);  // returned in kB
		pclose(fp2);
		if (ramkb > 0) {
			return ramkb * 1000.;
		}
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
			long long free_memory = ((int64_t)vm_stats.free_count +
                               (int64_t)vm_stats.inactive_count) * (int64_t)page_size;
			ram = free_memory;
		//https://stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process
		}
	#endif
	
	return ram;
}
