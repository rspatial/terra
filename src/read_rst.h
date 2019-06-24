std::vector<double> readBinAll(std::string filename, std::string datatype, 
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order);
							   
std::vector<double> readBinRows(std::string filename, std::string datatype, 
                               unsigned row, unsigned nrows,
                               unsigned nr, unsigned nc, unsigned nl,
                               std::string order);
					   

std::vector<double> readBinBlock(std::string filename, std::string datatype, 
                                unsigned row, unsigned nrows,
                                unsigned col, unsigned ncols,
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order);					
								
std::vector<double> readBinCell(std::string filename, std::string datatype, 
                                std::vector<double> cells,
                                unsigned nr, unsigned nc, unsigned nl,
                                std::string order);
