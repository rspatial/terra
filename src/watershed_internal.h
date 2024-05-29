// C/C++ code
// Author: Ezio Crestaz,Emanuele Cordano
// Date: February 2021
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer 1D array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function recursively call itself up to the watershed boundary or raster limit
// 
//


// Function: offset
// Scope: return offset of a raster cell stored in a 1D integers array with respect to base address
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
//
// not exported [[Rcpp::export]]
int offset(int nx, int ny, int x, int y);
//
// Functions: getRow, getCol
// Scope: return offset of a raster cell stored in a 1D integers array with respect to base address
// nx,ny: number of cells in x (longitude) and in y (latitude)
// offset: 
// NOTE: it assumes the offset is a valid value
//

int getRow(int nx, int ny, int offset);

int getCol(int nx, int ny, int offset);


//  Function: inRaster
// Scope: check if a given cell is valid or not, that's within the raster
//  nx,ny: number of cells in x (longitude) and in y (latitude)
//  x,y: indexes of cell upstream of which the watershed must be computed
// 

bool inRaster(int nx, int ny, int x, int y);
//  Function: watershed (version 0, recursive)

void watershed(double* p, int nx, int ny, int pp_offset, int* pOut);

//  Function: watershed (version 1, see below)
//
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Cells to be processed are dinamically managed in a queue, up to the basin boundaries
//
//void watershed_v1(double* p, int nx, int ny, int x, int y, int* pOut);
void watershed_v1(double* p, int nx, int ny, int x, int y, double* pOut);

//
// Function:   resizeQueue
// Scope:      resize (doubling) an existing queue to store integer offsets of raster cells
// Parameters: int* q pointer to existing queue
//             int n  current number of elements
// Return a pointer to the resized queue, while preserving previous values
//
int* resizeQueue(int* q, int n);

//
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Current version implements a dinamic queue, up to the basin boundaries
//
// void watershed_v2(int* p, int nx, int ny, int x, int y, int* pOut)
void watershed_v2(double* p, int nx, int ny, int pp_offset, double* pOut);



void pitfinder(double* p, int nx, int ny, double* pOut);
