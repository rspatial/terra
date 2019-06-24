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
                            unsigned nc, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(1);
  size_t n = cells.size();
  std::vector<T> out(n * nl);
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<n; i++) {
      size_t row = std::floor(cells[i] / nc);
      size_t col = cells[i] - (row * nc);
      size_t basestart = (row * nc * nl + col) * size;
      for (size_t j=0; j<nl; j++) {
        size_t start = basestart + (j * nc) * size;
        ifs.seekg(start, std::ios::beg);
        ifs.read(reinterpret_cast<char*>(d.data()), size);
        out[i + n*j] = d[0];
      }
  }
  return out;
}


template <typename T>
std::vector<T> readCellBSQ(std::string filename, 
                           std::vector<double> cells, 
                           unsigned nr, unsigned nc, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(1);
  size_t n = cells.size();
  size_t offset = nr * nc * size;
  std::vector<T> out(n * nl);
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<n; i++) {
    size_t row = std::floor(cells[i] / nc);
    size_t col = cells[i] - (row * nc);
    size_t basestart = (row * nc + col) * size;
    for (size_t j=0; j<nl; j++) {
      size_t start = basestart + offset;
      ifs.seekg(start, std::ios::beg);
      ifs.read(reinterpret_cast<char*>(d.data()), size);
      out[i + n*j] = d[0];
    }
  }
  return out;
}


template <typename T>
std::vector<T> readAll(std::string filename){
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  ifs.seekg(0, std::ios::end);
  std::streampos fileSize = ifs.tellg();
  ifs.seekg(0, std::ios::beg);
  size_t sizeOfBuffer = fileSize / sizeof(T);
  std::vector<T> d(sizeOfBuffer);
  ifs.read(reinterpret_cast<char*>(d.data()), fileSize);    
  return d;
}


template <typename T>
std::vector<T> readRowsBIL(std::string filename, unsigned row, unsigned nrows, unsigned nc, unsigned nl){
  size_t size = sizeof(T);
  size_t start = row * nc * nl * size;
  size_t n = nrows * nc * nl;
  
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  ifs.seekg(start, std::ios::beg);
  std::vector<T> d(n);
  ifs.read(reinterpret_cast<char*>(d.data()), n * size);    
  return d;
}


template <typename T>
std::vector<T> readRowsBSQ(std::string filename, unsigned row, unsigned nrows, unsigned nr, unsigned nc, unsigned nl){
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
  return out;
}

template <typename T>
std::vector<T> readBlockBIL(std::string filename, 
                            unsigned row, unsigned nrows, 
                            unsigned col, unsigned ncols, 
                            unsigned nc, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(ncols);
  std::vector<T> out;
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<nl; i++) {
    for (size_t j=0; j<nrows; j++) {
      size_t start = ((row+j) * nc * nl + i * nc + col) * size;
      ifs.seekg(start, std::ios::beg);
      ifs.read(reinterpret_cast<char*>(d.data()), ncols * size);
      out.insert(out.end(), d.begin(), d.end());
    }
  }
  return out;
}


template <typename T>
std::vector<T> readBlockBSQ(std::string filename, 
                            unsigned row, unsigned nrows, 
                            unsigned col, unsigned ncols, 
                            unsigned nc, unsigned nr, unsigned nl){
  size_t size = sizeof(T);
  std::vector<T> d(ncols);
  std::vector<T> out;
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  for (size_t i=0; i<nl; i++) {
    for (size_t j=0; j<nrows; j++) {
      size_t start = (nr * nc * i + (row+j) * nc + col) * size;
      ifs.seekg(start, std::ios::beg);
      ifs.read(reinterpret_cast<char*>(d.data()), ncols * size);
      out.insert(out.end(), d.begin(), d.end());
    }
  }
  return out;
}



std::vector<double> readBinRows(std::string filename, std::string datatype, 
                               unsigned row, unsigned nrows,
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order) {
  std::vector<double> out;
  if (order == "BIL") {

    if (datatype == "INT2S") {
      std::vector<short> v = readRowsBIL<short>(filename, row, nrows, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readRowsBIL<long>(filename, row, nrows, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readRowsBIL<float>(filename, row, nrows, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readRowsBIL<double>(filename, row, nrows, nc, nl);  
    }
    out = bil_to_bsq(out, nrows, nc, nl );
  } else if (order == "BSQ") {
    if (datatype == "INT2S") {
      std::vector<short> v = readRowsBSQ<short>(filename, row, nrows, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readRowsBSQ<long>(filename, row, nrows, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readRowsBSQ<float>(filename, row, nrows, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readRowsBSQ<double>(filename, row, nrows, nr, nc, nl);  
    }
  }
  return out;
}


std::vector<double> readBinBlock(std::string filename, std::string datatype, 
                                unsigned row, unsigned nrows,
                                unsigned col, unsigned ncols,
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order) {
  std::vector<double> out;
  if (order == "BIL") {
    
    if (datatype == "INT2S") {
      std::vector<short> v = readBlockBIL<short>(filename, row, nrows, col, ncols, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readBlockBIL<long>(filename, row, nrows, col, ncols, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readBlockBIL<float>(filename, row, nrows, col, ncols, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readBlockBIL<double>(filename, row, nrows, col, ncols, nc, nl);  
    }
    out = bil_to_bsq(out, nrows, ncols, nl );
  } else if (order == "BSQ") {
    if (datatype == "INT2S") {
      std::vector<short> v = readBlockBSQ<short>(filename, row, nrows, col, ncols, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readBlockBSQ<long>(filename, row, nrows, col, ncols, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readBlockBSQ<float>(filename, row, nrows, col, ncols, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readBlockBSQ<double>(filename, row, nrows, col, ncols, nr, nc, nl);  
    }
  }
  return out;
}


std::vector<double> readBinAll(std::string filename, std::string datatype, 
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order) {
    std::vector<double> out;
    if (datatype == "INT2S") {
      std::vector<short> v = readAll<short>(filename); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readAll<long>(filename); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readAll<float>(filename); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readAll<double>(filename);  
    }
    if (order == "BIL") {
      out = bil_to_bsq(out, nr, nc, nl );
    }
    return out;
}



std::vector<double> readBinCell(std::string filename, std::string datatype, 
                                std::vector<double> cells,
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order) {
  std::vector<double> out;
  if (order == "BIL") {
    
    if (datatype == "INT2S") {
      std::vector<short> v = readCellBIL<short>(filename, cells, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readCellBIL<long>(filename, cells, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readCellBIL<float>(filename, cells, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readCellBIL<double>(filename, cells, nc, nl);  
    }
  } else if (order == "BSQ") {
    if (datatype == "INT2S") {
      std::vector<short> v = readCellBSQ<short>(filename, cells, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "INT4S") {
      std::vector<long> v = readCellBSQ<long>(filename, cells, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT4S") {
      std::vector<float> v = readCellBSQ<float>(filename, cells, nr, nc, nl); 
      out = std::vector<double>(v.begin(), v.end());
    } else if (datatype == "FLT8S") {
      out = readCellBSQ<double>(filename, cells, nr, nc, nl);  
    }
  }
  return out;
}
