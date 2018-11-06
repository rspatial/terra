#include <vector>
#include <string>

class SpatDataFrame {
	public:
		std::vector<std::string> names;
		std::vector<unsigned> itype;
		std::vector<unsigned> iplace;
		std::vector< std::vector<double>> dv;
		std::vector< std::vector<long>> iv;
		std::vector< std::vector<std::string>> sv;
		std::string NAS = "____NA_+";
		
		unsigned nrow();
		unsigned ncol();
		SpatDataFrame subset_rows(std::vector<unsigned> range);
		SpatDataFrame subset_cols(std::vector<unsigned> range);
		SpatDataFrame subset_rows(unsigned i);
		SpatDataFrame subset_cols(unsigned i);
		std::vector<double> getD(unsigned i);
		std::vector<long> getI(unsigned i);
		std::vector<std::string> getS(unsigned i);
	
		void add_row();
		void add_column(unsigned dtype, std::string name);
		bool add_column(std::vector<double> x, std::string name);
		bool add_column(std::vector<long> x, std::string name);
		bool add_column(std::vector<std::string> x, std::string name);
		
		void insert_column(std::vector<double>);
		void insert_column(std::vector<long>);		
		void insert_column(std::vector<std::string>);

		void resize(unsigned n);
		void reserve(unsigned n);
		
		bool rbind(SpatDataFrame &x);
		bool cbind(SpatDataFrame &x);
		
};

