using namespace std;
#include <vector>
#include <string>

class SpatDataFrame {
	public:
		std::vector<string> names;
		std::vector<unsigned> itype;
		std::vector<unsigned> iplace;
		std::vector< std::vector<double>> dv;
		std::vector< std::vector<long>> iv;
		std::vector< std::vector<string>> sv;
		string NAS = "____NA_+";
		
		unsigned nrow();
		unsigned ncol();
		SpatDataFrame subsetrows(std::vector<unsigned> range);
		SpatDataFrame subsetcols(std::vector<unsigned> range);
		std::vector<double> getD(unsigned i);
		std::vector<long> getI(unsigned i);
		std::vector<string> getS(unsigned i);
	
	void add_row();
};

