#include "spatRaster.h"

std::vector<std::string> unpack(unsigned value, int nbits, std::vector<unsigned> idx) {
  
  std::vector<int> b(nbits); 
  for (int i = 0; i < nbits; i++) {
    b[nbits-i-1] = (value >> i) & 1;
  }
  
  std::vector<std::string> s;
  for (size_t i = 0; i < idx.size(); i=i+2) {
    if (idx[i] == idx[i+1]) {
      s.push_back(std::to_string(b[idx[i]]));
    } else {
      std::vector<int> slice(b.begin()+idx[i], b.begin()+idx[i+1]+1);
      if (slice.size() > 1) {
        std::string ss = "";
        for(int& j : slice)  ss = ss + std::to_string(j);
        s.push_back(ss);
      } 
    }
  }  
  return s;  
} 
  
std::vector<int> matchbits(std::vector<unsigned> values, int nbits, std::vector<unsigned> idx, std::vector<std::string> match, std::unordered_map<unsigned, bool> &umap) {
	std::vector<int> matches(values.size());
	for (size_t i=0; i<values.size(); i++) {
		unsigned v = values[i];
		if (umap.find(v) == umap.end()) {
			std::vector<std::string> s = unpack(v, nbits, idx);
			for (size_t j=0; j<s.size(); j++) {
				umap[v] = 0;
				if (s[j] == match[j]) {
					matches[i] = 1; 
					umap[v] = 1;
					break;
				}
			}	
		} else {
			matches[i] = umap[v];
		}
	}  
	return matches; 
}



SpatRaster SpatRaster::modisqc(int nbits, std::vector<unsigned> idx, std::vector<std::string> match, SpatOptions &opt) {
	SpatRaster out = geometry();
	if (!hasValues()) {
		out.setError("qc layer has no values");
		return out;
	}
	if (nlyr() > 1) {
		out.setError("qc object has multiple layers");
		return out;		
	}
	
	std::unordered_map<unsigned, bool> umap; 
	readStart();
 	if (!out.writeStart(opt)) { return out; }
	for (size_t i = 0; i < out.bs.n; i++) {
		std::vector<double> v = readBlock(out.bs, i);
		std::vector<unsigned> vv(v.begin(), v.end());
		std::vector<int> m = matchbits(vv, nbits, idx, match, umap);
		std::vector<double> mm(m.begin(), m.end());
		if (!out.writeValues(mm, i)) return out;
	}
	readStop();
	out.writeStop();
	return(out);
}

