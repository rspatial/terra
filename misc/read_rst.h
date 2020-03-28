
std::vector<double> readBinAll(std::string filename, std::string datatype,
							   std::vector<unsigned> lyrs,
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order);
							   
std::vector<double> readBinRows(std::string filename, std::string datatype, 
                               unsigned row, unsigned nrows,
							   std::vector<unsigned> lyrs,
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order);
					   
std::vector<double> readBinBlock(std::string filename, std::string datatype, 
                                unsigned row, unsigned nrows,
                                unsigned col, unsigned ncols,
 							    std::vector<unsigned> lyrs,
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order);					
								
std::vector<std::vector<double>> readBinCell(std::string filename, std::string datatype, 
                                std::vector<double> cells,
							    std::vector<unsigned> lyrs,
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order);

