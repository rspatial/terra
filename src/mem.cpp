#include "spat.h"

#ifdef _WIN32
#include <windows.h>
#else
#include "sys/types.h"
#include "sys/sysinfo.h"
#include <unistd.h>
#endif


double availableRAM() {
	// return available RAM in number of double (8 byte) cells.
	double ram;
	#ifdef _WIN32
		MEMORYSTATUSEX statex;
		statex.dwLength = sizeof(statex);
		GlobalMemoryStatusEx(&statex);
		ram = statex.ullAvailPhys / 8.;
	#else
		struct sysinfo memInfo;
		sysinfo (&memInfo);
		ram = memInfo.freeram / 8.;
	#endif
	
	return ram;
}


bool SpatRaster::canProcessInMemory(unsigned n) {
	double f = 0.5;
	return (n * size()) < (availableRAM() * f);
}



unsigned SpatRaster::chunkSize(unsigned n) {
	double f = 0.25;
	unsigned cells_in_row = n * ncol * nlyr();
	unsigned rows = availableRAM() * f / cells_in_row;
	return rows == 0 ? 1 : min(rows, nrow);
}

BlockSize SpatRaster::getBlockSize(unsigned n) {
	BlockSize bs;

//	if (source[0].filename == "") {
	// in memory
//		bs.row = {0};
//		bs.nrows = {nrow};
//		bs.n = 1;

//	} else {

		unsigned cs = chunkSize(n);
		unsigned chunks = ceil(nrow / double(cs));
		bs.n = chunks;
		bs.row = vector<unsigned>(chunks, 0);
		bs.nrows = vector<unsigned>(chunks, cs);

		unsigned r = 0;
		for (size_t i =0; i<chunks; i++) {
			bs.row[i] = r;
			r += cs;
		}	
		bs.nrows[chunks-1] = cs - ((chunks * cs) - nrow);	

//	}
	return bs;
}
