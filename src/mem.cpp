#include "spat.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef linux
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
	double f = 0.5;
	unsigned rows = availableRAM() * f / n / ncol;
	return min(rows, nrow);
}

BlockSize SpatRaster::getBlockSize() {
	BlockSize bs;

	if (source[0].filename == "") {
	// in memory
		bs.row = {0};
		bs.nrows = {nrow};
		bs.n = 1;

	} else {

//		unsigned rowsperchunk = chunkSize();

	// to be improved
		bs.row = {0, unsigned(floor(nrow/2))};
		bs.nrows = {bs.row[1], nrow-bs.row[1]};
		bs.n = 2;
	}
	return bs;
}
