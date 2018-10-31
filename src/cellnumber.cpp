using namespace std;

#include "spatraster.h"
#include <cmath>



std::vector<double> SpatRaster::cellFromXY (std::vector<double> x, std::vector<double> y) {

// size of x and y should be the same
	size_t size = x.size();
	std::vector<double> cells(size);
	
	double yr_inv = nrow / (extent.ymax - extent.ymin);
	double xr_inv = ncol / (extent.xmax - extent.xmin);
  
	for (size_t i = 0; i < size; i++) {
		// cannot use trunc here because trunc(-0.1) == 0
		unsigned row = std::floor((extent.ymax - y[i]) * yr_inv);
		// points in between rows go to the row below
		// except for the last row, when they must go up
		if (y[i] == extent.ymin) {  
			row = nrow-1 ;
		}
		
		unsigned col = floor((x[i] - extent.xmin) * xr_inv);
		// as for rows above. Go right, except for last column
		if (x[i] == extent.xmax) {
			col = ncol - 1 ;
		}
		
		if (row < 0 || row >= nrow || col < 0 || col >= ncol) {
			cells[i] = NAN;
		} else {
			// result[i] = static_cast<int>(row) * ncols + static_cast<int>(col) + 1;
			cells[i] = row * ncol + col;
		}
	}
  
	return cells;
}


double SpatRaster::cellFromXY (double x, double y) {
	std::vector<double> X = {x};
	std::vector<double> Y = {y};
	std::vector<double> cell = cellFromXY(X, Y);
	return  cell[0]; 
}


std::vector<double> SpatRaster::cellFromRowCol(std::vector<unsigned> rownr, std::vector<unsigned> colnr) {
	
	size_t rownr_size = rownr.size();
	size_t colnr_size = colnr.size();
  
	std::vector<double> result(std::max(rownr_size, colnr_size));

  // Manually recycle the shorter of rownr/colnr to match the other
	size_t len = std::max(rownr.size(), colnr.size());

	for (size_t i = 0; i < len; i++) {
    // The % is to recycle elements if they're not the same length
		double r = rownr[i < rownr_size ? i : i % rownr_size];
		double c = colnr[i < colnr_size ? i : i % colnr_size];

    // Detect out-of-bounds rows/cols and use NA for those
		result[i] = (r < 0 || r >= nrow || c < 0 || c >= ncol) ? NAN : r * ncol + c;
	}
  
	return result;
}


double SpatRaster::cellFromRowCol (unsigned rownr, unsigned colnr) {
	std::vector<unsigned> rows = {rownr};
	std::vector<unsigned> cols = {colnr};
	std::vector<double> cell = cellFromRowCol(rows, cols);
	return  cell[0]; 
}


std::vector<double> SpatRaster::yFromRow(std::vector<unsigned> rownr) {
	size_t size = rownr.size();
	std::vector<double> result( size );
	double ymax = extent.ymax;
	double yr = yres();
	for (size_t i = 0; i < size; i++) {
		result[i] = (rownr[i] < 0 || rownr[i] >= nrow ) ? NAN : ymax - ((rownr[i]+0.5) * yr);
	}
	return result;
}
	
double SpatRaster::yFromRow (unsigned rownr) {
	std::vector<unsigned> rows = {rownr};
	std::vector<double> y = yFromRow(rows);
	return y[0]; 
}



std::vector<double> SpatRaster::xFromCol(std::vector<unsigned> colnr) {
	size_t size = colnr.size();
	std::vector<double> result( size );
	double xmin = extent.xmin;
	double xr = xres();
	for (size_t i = 0; i < size; i++) {
		result[i] = (colnr[i] < 0 || colnr[i] >= ncol ) ? NAN : xmin + ((colnr[i]+0.5) * xr);
	}
	return result;
}
	
double SpatRaster::xFromCol(unsigned colnr) {
	std::vector<unsigned> cols = {colnr};
	std::vector<double> x = xFromCol(cols);
	return x[0]; 
}

std::vector<double> SpatRaster::colFromX(std::vector<double> x) {
	size_t size = x.size();
	std::vector<double> result(size);
	double xmin = extent.xmin;
	double xmax = extent.xmax;
	double xr = xres();
	
	for (size_t i = 0; i < size; i++) {
		if (x[i] == xmax) {  
			result[i] = ncol ;
		} else {
			result[i] = (x[i] < xmin || x[i] > xmax ) ? NAN : trunc((x[i] - xmin) / xr);
		}
	}
	return result;
}


double SpatRaster::colFromX(double x) {
	std::vector<double> X = {x};
	return colFromX(X)[0];
}


std::vector<double> SpatRaster::rowFromY(std::vector<double> y) {
	size_t size = y.size();
	std::vector<double> result(size);
	double ymin = extent.ymin;
	double ymax = extent.ymax;
	double yr = yres();
	
	for (size_t i = 0; i < size; i++) {
		if (y[i] == ymin) {  
			result[i] = nrow ;
		} else {
			result[i] = (y[i] < ymin || y[i] > ymax ) ? NAN : trunc((ymax - y[i]) / yr);
		}
	}
	return result;
}

double SpatRaster::rowFromY(double y) {
	std::vector<double> Y = {y};
	return rowFromY(Y)[0];
}



std::vector< std::vector<double> > SpatRaster::xyFromCell( std::vector<double> cell ) {
	size_t size = cell.size();
	double xmin = extent.xmin;
	double ymax = extent.ymax;
	double yr = yres();
	double xr = xres();
  
	std::vector< std::vector<double> > result(2, std::vector<double> (size) );
	for (size_t i = 0; i < size; i++) {
		unsigned row = (cell[i] / ncol);
		unsigned col = fmod(cell[i], ncol);
		//unsigned col = cell[i] - row * ncol
		result[0][i] = xmin + (col + 0.5) * xr;
		result[1][i] = ymax - (row + 0.5) * yr;
	}
	return result;
}


std::vector< std::vector<double> > SpatRaster::xyFromCell( double cell ) {
	std::vector<double> Cell = {cell};
	return xyFromCell(Cell);
}


std::vector< std::vector<double> > SpatRaster::rowColFromCell(std::vector<double> cell) {
	size_t size = cell.size();
	std::vector< std::vector<double> > result(2, std::vector<double> (size) );

	double nc = ncell();
	
	for (size_t i = 0; i < size; i++) {
		if ((cell[i] < 0 || cell[i] >= nc )) {  
			result[0][i] = NAN;
			result[1][i] = NAN;
		} else {
			result[0][i] = trunc(cell[i]/ ncol);
			result[1][i] = (cell[i] - ((result[0][i]) * ncol));
		}
	}
	return result;
}


