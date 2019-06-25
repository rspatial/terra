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
#include <cmath> // floor

/*
https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
#include <climits>
template <typename T>
T swap_endian(T u) {
    static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");
    union {
        T u;
        unsigned char u8[sizeof(T)];
    } source, dest;

    source.u = u;
    for (size_t k = 0; k < sizeof(T); k++) {
        dest.u8[k] = source.u8[sizeof(T) - k - 1];
    }
    return dest.u;
}
swap_endian<double>(42).
*/

template <typename T>
std::vector<T> bil_to_bsq(const std::vector<T> &v, unsigned nrows, unsigned ncols, unsigned nlyrs) {
  std::vector<T> x;
  for (size_t i=0; i<nlyrs; i++) {
    for (size_t r=0; r<nrows; r++) {
      unsigned step = (r * nlyrs + i) * ncols;
      x.insert(x.end(), v.begin()+step, v.begin()+step+ncols);
    }
  }
  return x;
}


template <typename T>
std::vector<T> readCellBIL(std::string filename, 
                           std::vector<double> cells, 
                           std::vector<unsigned> lyrs,
                           unsigned nc, unsigned nl){
  
	size_t size = sizeof(T);
	std::vector<T> d(1);
	size_t n = cells.size();
	size_t nlyrs = lyrs.size();
	std::vector<T> out(n * nlyrs);
	std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  
	for (size_t i=0; i<n; i++) {
		size_t row = std::floor(cells[i] / nc);
		size_t col = cells[i] - (row * nc);
		size_t basepos = (row * nc * nl + col) * size;
		for (size_t j=0; j<nlyrs; j++) {
			size_t pos = basepos + lyrs[j] * nc * size;
			ifs.seekg(pos, std::ios::beg);
			ifs.read(reinterpret_cast<char*>(d.data()), size);
			out[i + n*j] = d[0];
		}
	}
	ifs.close();
	return out;
}


template <typename T>
std::vector<T> readCellBSQ(std::string filename, 
                           std::vector<double> cells, 
                           std::vector<unsigned> lyrs,
                           unsigned nr, unsigned nc, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(1);
  size_t n = cells.size();
  size_t nlyrs = lyrs.size();
  size_t offset = nr * nc * size;
  std::vector<T> out(n * nlyrs);
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<n; i++) {
    size_t row = std::floor(cells[i] / nc);
    size_t col = cells[i] - (row * nc);
    size_t basepos = (row * nc + col) * size;
    for (size_t j=0; j<nlyrs; j++) {
      size_t pos = basepos + offset * j;
      ifs.seekg(pos, std::ios::beg);
      ifs.read(reinterpret_cast<char*>(d.data()), size);
      out[i + n*j] = d[0];
    }
  }
  ifs.close();
  return out;
}

template <typename T>
std::vector<T> readAll(std::string filename, std::vector<unsigned> lyrs, unsigned nr, unsigned nc, unsigned nl, std::string bandorder){
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  ifs.seekg(0, std::ios::end);
  std::streampos fileSize = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  size_t n = fileSize / sizeof(T);
  std::vector<T> d(n);
  ifs.read(reinterpret_cast<char*>(d.data()), fileSize);    
  ifs.close();

  if (bandorder == "BIL") {
    d = bil_to_bsq(d, nr, nc, nl);
  }
  
  size_t nlyrs = lyrs.size();
  if (nlyrs != nl) {
    size_t ncell = nc * nr;
    std::vector<T> out;
    for (size_t i=0; i<nlyrs; i++) {
      size_t start = lyrs[i] * ncell;
      out.insert(out.end(), d.begin()+start, d.begin()+start+ncell);
    }	
    d = out;
  }
  return d;
}


template <typename T>
std::vector<T> readRowsBIL(std::string filename, 
                           unsigned row, unsigned nrows, 
                           std::vector<unsigned> lyrs, 
                           unsigned nc, unsigned nl) {
  
  size_t size = sizeof(T);
  size_t start = row * nc * nl * size;
  size_t n = nrows * nc * nl;
  
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  ifs.seekg(start, std::ios::beg);
  std::vector<T> d(n);
  ifs.read(reinterpret_cast<char*>(d.data()), n * size);    
  ifs.close();

  d = bil_to_bsq(d, nrows, nc, nl);
  
  size_t nlyrs = lyrs.size();
  if (nlyrs != nl) {
    std::vector<T> out;
    for (size_t i=0; i<nlyrs; i++) {
      size_t start = lyrs[i] * nrows * nc;
      out.insert(out.end(), d.begin()+start, d.begin()+start+nrows*nc);
    }	
    d = out;
  }
  
  return d;
}


template <typename T>
std::vector<T> readRowsBSQ(std::string filename, 
                           unsigned row, unsigned nrows, 
                           std::vector<unsigned> lyrs, 
                           unsigned nr, unsigned nc, unsigned nl){
  
  size_t size = sizeof(T);
  size_t n = nrows * nc;
  std::vector<T> d(n);
  std::vector<T> out;
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<nl; i++) {
    size_t start = (row * nc + (i * nr)) * size;
    ifs.seekg(start, std::ios::beg);
    ifs.read(reinterpret_cast<char*>(d.data()), n * sizeof(T));   
    out.insert(out.end(), d.begin(), d.end());
  }
  ifs.close();
  return out;
}

template <typename T>
std::vector<T> readBlockBIL(std::string filename, 
                            unsigned row, unsigned nrows, 
                            unsigned col, unsigned ncols, 
                            std::vector<unsigned> lyrs, 
                            unsigned nc, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(ncols);
  std::vector<T> out;
  size_t nlyrs = lyrs.size();
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<nlyrs; i++) {
    for (size_t j=0; j<nrows; j++) {
      size_t start = ((row+j) * nc * nl + lyrs[i] * nc + col) * size;
      ifs.seekg(start, std::ios::beg);
      ifs.read(reinterpret_cast<char*>(d.data()), ncols * size);
      out.insert(out.end(), d.begin(), d.end());
    }
  }
  ifs.close();
  out = bil_to_bsq(out, nrows, ncols, nlyrs);
  return out;
}


template <typename T>
std::vector<T> readBlockBSQ(std::string filename, 
                            unsigned row, unsigned nrows, 
                            unsigned col, unsigned ncols, 
                            std::vector<unsigned> lyrs, 
                            unsigned nc, unsigned nr, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(ncols);
  std::vector<T> out;
  size_t nlyrs = lyrs.size();
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<nlyrs; i++) {
    for (size_t j=0; j<nrows; j++) {
      size_t start = (nr * nc * lyrs[i] + (row+j) * nc + col) * size;
      ifs.seekg(start, std::ios::beg);
      ifs.read(reinterpret_cast<char*>(d.data()), ncols * size);
      out.insert(out.end(), d.begin(), d.end());
    }
  }
  ifs.close();
  return out;
}



std::vector<double> readBinRows(std::string filename, std::string datatype, 
                                unsigned row, unsigned nrows,
                                std::vector<unsigned> lyrs, 
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order) {
  std::vector<double> out;
  if (order == "BIL") {
    
    if (datatype == "INT2S") {
      std::vector<short> v = readRowsBIL<short>(filename, row, nrows, lyrs, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readRowsBIL<long>(filename, row, nrows, lyrs, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readRowsBIL<float>(filename, row, nrows, lyrs, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readRowsBIL<double>(filename, row, nrows, lyrs, nc, nl);  
    }
  } else if (order == "BSQ") {
    if (datatype == "INT2S") {
      std::vector<short> v = readRowsBSQ<short>(filename, row, nrows, lyrs, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readRowsBSQ<long>(filename, row, nrows, lyrs, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readRowsBSQ<float>(filename, row, nrows, lyrs, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readRowsBSQ<double>(filename, row, nrows, lyrs, nr, nc, nl);  
    }
  }
  return out;
}


std::vector<double> readBinBlock(std::string filename, std::string datatype, 
                                 unsigned row, unsigned nrows,
                                 unsigned col, unsigned ncols,
                                 std::vector<unsigned> lyrs, 
                                 unsigned nr, unsigned nc, unsigned nl,
                                 std::string order) {
  std::vector<double> out;
  if (order == "BIL") {
    
    if (datatype == "INT2S") {
      std::vector<short> v = readBlockBIL<short>(filename, row, nrows, col, ncols, lyrs, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readBlockBIL<long>(filename, row, nrows, col, ncols, lyrs, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readBlockBIL<float>(filename, row, nrows, col, ncols, lyrs, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readBlockBIL<double>(filename, row, nrows, col, ncols, lyrs, nc, nl);  
    }
  } else if (order == "BSQ") {
    if (datatype == "INT2S") {
      std::vector<short> v = readBlockBSQ<short>(filename, row, nrows, col, ncols, lyrs, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readBlockBSQ<long>(filename, row, nrows, col, ncols, lyrs, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readBlockBSQ<float>(filename, row, nrows, col, ncols, lyrs, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readBlockBSQ<double>(filename, row, nrows, col, ncols, lyrs, nr, nc, nl);  
    }
  }
  return out;
}


std::vector<double> readBinAll(std::string filename, std::string datatype, 
                               std::vector<unsigned> lyrs, 
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order) {
  std::vector<double> out;
  if (datatype == "INT2S") {
    std::vector<short> v = readAll<short>(filename, lyrs, nr, nc, nl, order); 
    out = std::vector<double>(v.begin(), v.end());
  } else if (datatype == "INT4S") {
    std::vector<long> v = readAll<long>(filename, lyrs, nr, nc, nl, order); 
    out = std::vector<double>(v.begin(), v.end());
  } else if (datatype == "FLT4S") {
    std::vector<float> v = readAll<float>(filename, lyrs, nr, nc, nl, order); 
    out = std::vector<double>(v.begin(), v.end());
  } else if (datatype == "FLT8S") {
    out = readAll<double>(filename, lyrs, nr, nc, nl, order);  
  }
  return out;
}


std::vector<std::vector<double>> readBinCell(std::string filename, 
			std::string datatype, std::vector<double> cells,
            std::vector<unsigned> lyrs,
            unsigned nr, unsigned nc, unsigned nl, std::string order) {
				
	std::vector<double> d;
	
	if (order == "BIL") {    
		if (datatype == "INT2S") {
			std::vector<short> v = readCellBIL<short>(filename, cells, lyrs, nc, nl); 
			d = std::vector<double>(v.begin(), v.end());
		} else if (datatype == "INT4S") {
			std::vector<long> v = readCellBIL<long>(filename, cells, lyrs, nc, nl); 
			d = std::vector<double>(v.begin(), v.end());
		} else if (datatype == "FLT4S") {
			std::vector<float> v = readCellBIL<float>(filename, cells, lyrs, nc, nl); 
			d = std::vector<double>(v.begin(), v.end());
		} else if (datatype == "FLT8S") {
			d = readCellBIL<double>(filename, cells, lyrs, nc, nl);  
		}
	} else if (order == "BSQ") {
		if (datatype == "INT2S") {
			std::vector<short> v = readCellBSQ<short>(filename, cells, lyrs, nr, nc, nl); 
			d = std::vector<double>(v.begin(), v.end());
		} else if (datatype == "INT4S") {
			std::vector<long> v = readCellBSQ<long>(filename, cells, lyrs, nr, nc, nl); 
			d = std::vector<double>(v.begin(), v.end());
		} else if (datatype == "FLT4S") {
			std::vector<float> v = readCellBSQ<float>(filename, cells, lyrs, nr, nc, nl); 
			d = std::vector<double>(v.begin(), v.end());
		} else if (datatype == "FLT8S") {
			d = readCellBSQ<double>(filename, cells, lyrs, nr, nc, nl);  
		}
	}
  
	size_t nlyrs = lyrs.size();
	size_t n = cells.size();

	std::vector<std::vector<double>> out(nlyrs, std::vector<double>(n)) ;
	for (size_t i=0; i<nlyrs; i++) {
		size_t start = i * n;
		out[i] = std::vector<double>(d.begin()+start, d.begin()+start+n);
	}
    return out;
}


