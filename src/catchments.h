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
void pitfinder(double* p, int nx, int ny, double* pOut,int pits_on_boundary);

///// flow accumulation
void NextCell(double* p, int nx, int ny,int* pnext);
void NIDP(int* pnext, int nx, int ny,double* nidp_value); 
// FlowAccu algorithm 5 
// Reference: https://link.springer.com/article/10.1007/s11707-018-0725-9
void FlowAccu(int* pnext, int nx, int ny,double* nidp_value,double* flowaccu_value);

void FlowAccu_weight(int* pnext, int nx, int ny,double* nidp_value,double* flowaccu_value,double* weight);

int nextcell_point_conv1(int nx, int ny,int x,int y,double pv,int conv_type);


// 
// void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,
//                      double* tdc, double *tdd,double L,
//                      std::vector<double> ddp1,std::vector<double> ddp2,int nncell,int conv_type,int use_lad);

void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,
                     double* tdc, double *tdd,double L,
                     std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,
                     int conv_type,int use_lad);


void transverse_deviation(double *e, double *tdc, double *tdd,double *sr,double *sm, int *sfacet,int nx, int ny, double L,
                          
                          double *atdc, double *atdd, double *atdplus, double *atdplus0,
                          double *pflow,int *has_upstream,int *kupdate,double *npids,
                          double lambda,
                          std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,int conv_type,int use_lad,int max_iters);

void d8ltd_computation(double *e,int nx,int ny,double L,double lambda,int use_lad,int max_iters,double *pflow);




//// pitfiller 
#define NAN_PITFILLER -9999.9
#define WINIT_JK -1
#define FALSE_FILLED 0
#define TRUE_FILLED 1

////////void pitfiller(int nx,int ny,double *pitf,double* e,double *eout,double*flowdir,double L);
void pitfiller(int nx,int ny,double ipit, double *pitf,double *e,double *eout,double *flowdir,double L,double *flowacc,
               double U,double D,double beta,double theta_exp); // see reference doi:10.1016/j.advwatres.2006.11.016
/////void pitfiller_all(int nx,int ny, double *pitf,double *pitfdone,double* e,double *eout,double *flowdir); 

void pitfiller_all(int nx,int ny, double *pitf,double *pitftemp,double* e,double *eout,double *flowdir,int niter,double lambda,int use_lad,double L,
                    double U,double D,double beta,double theta_exp); // see // see reference doi:10.1016/j.advwatres.2006.11.016)    


double calculeZ_pem4ppt(double U, double delta_l, double dx, double D, double z_avg, double beta, double A, double theta, double z_d); 
