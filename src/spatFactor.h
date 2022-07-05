#ifndef SPATFACTOR_GUARD
#define SPATFACTOR_GUARD


class SpatFactor {
public:
	virtual ~SpatFactor(){}
	SpatFactor(){} ;
	SpatFactor(size_t _size, unsigned _value) {
		v.resize(_size, _value);
	};

	SpatFactor(std::vector<unsigned> _values, std::vector<std::string> _labels);
	SpatFactor(std::vector<unsigned> _values);
	SpatFactor(std::vector<std::string> _values);

	
	std::vector<unsigned> v;
	//std::vector<unsigned> levels;
	std::vector<std::string> labels;
  
	size_t size() { return v.size(); }
	
	//void compute_levels();
	void push_back(unsigned x) { v.push_back(x); }

	
	bool set_labels(std::vector<std::string> _labels);
	
	void reserve(size_t n) { v.reserve(n); }
	void resize(size_t n) { v.resize(n); }
	void resize(size_t n, unsigned x) {v.resize(n, x);}	
	
//	template <typename T>
//	  SpatFactor(std::vector<T> _v) {
//	   set_values(_v);
//	}
	
	SpatFactor subset(std::vector<unsigned> i);
	std::string getLabel(size_t i); 
	std::vector<std::string> getLabels();
	
};

#endif

