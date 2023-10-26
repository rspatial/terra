#include <Rcpp.h>
using namespace Rcpp;
#include "watershed_internal.h"
#include "spatRaster.h"




// C/C++ code
// Author: Ezio Crestaz
// Data: Febrary 2021
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
int offset(int nx, int ny, int x, int y)
{
  return y * nx + x; //according to original Ezio's code BY ROWS 
  //return x * ny + y; // according offset defintion in Rccp for IntergerMatrix // BY COLS
}



//
// Functions: getRow, getCol
// Scope: return offset of a raster cell stored in a 1D integers array with respect to base address
// nx,ny: number of cells in x (longitude) and in y (latitude)
// offset: 
// NOTE: it assumes the offset is a valid value
//



int getRow(int nx, int ny, int offset)

{
  return offset/nx;// according to Ezio's original code // BY ROWS
  //return offset % ny; // according offset definition in Rccp for IntergerMatrix // BY COLS
}

int getCol(int nx, int ny, int offset)// according to Ezio's original code

{
  return offset % nx;// according to Ezio's original code // BY ROWS
  //return offset/ny; // according offset definition in Rccp for IntergerMatrix // BY COLS
}










//  Function: inRaster
// Scope: check if a given cell is valid or not, that's within the raster
//  nx,ny: number of cells in x (longitude) and in y (latitude)
//  x,y: indexes of cell upstream of which the watershed must be computed
// 

bool inRaster(int nx, int ny, int x, int y)
{
  if (x < 0 || x >= nx || y >= ny || y < 0) return false;
  return true;
}


//  Function: watershed (version 0, recursive)

void watershed(double* p, int nx, int ny, int x, int y, int* pOut)
{
  static int nCall = 0;
  nCall++;
  ///printf("%d \n", nCall);
  ///printf("%d %d %d %d \n",nx,ny,x,y);
  
  // Set cell under analysis to 1, being part of the watershed
  // For the first cell it is not granted that it is a valid raster cell
  *(pOut + offset(nx, ny, x, y)) = 1;
  
  // Given the cell under analysis at x,y all 8 bounding cells are processed, checking the following: 
  // 1. the bounding cell must be within the raster, which is not true when the cell under analysis is 
  //    on the boundary of the input raster);
  // 2. the bounding cell in the output raster, that will identify the watershed upstream of the pour point, 
  //    must be 0, otherwise it means that it has already been identified as belonging to the watershed and
  //    hence it does not require any further processing;
  // 3. the bounding cell from the flow direction raster must report a flow direction towards the cell under
  //    analysis. The value change consistently with expected flow direction and ESRI ArcGIS codification
  // If, and only if, all the above conditions hold true, then the function watershed is called recursively
  // on the bounding cell. The watershed function itself will return back to the calling function when all
  // the upstream tree will be investigated and the watershed fully identified.
  
  // Bounding raster cell located to the E
  if (inRaster(nx, ny, x + 1, y) &&
  //    !*(pOut + offset(nx, ny, x + 1, y)) && 
  *(p + offset(nx, ny, x + 1, y)) == 16) 
    watershed(p, nx, ny, x + 1, y,pOut);
  // Bounding raster cell located to the SE
  if (inRaster(nx, ny, x + 1, y + 1) &&
      //    !*(pOut + offset(nx, ny, x + 1, y + 1)) && 
      *(p + offset(nx, ny, x + 1, y + 1)) == 32) 
    watershed(p, nx, ny, x + 1, y + 1, pOut);
  // Bounding raster cell located to the S
  if (inRaster(nx, ny, x, y + 1) &&
      //    !*(pOut + offset(nx, ny, x, y + 1)) && 
      *(p + offset(nx, ny, x, y + 1)) == 64) 
    watershed(p, nx, ny, x, y + 1, pOut);
  // Bounding raster cell located to the SW
  if (inRaster(nx, ny, x - 1, y + 1) && 
      //    !*(pOut + offset(nx, ny, x - 1, y + 1)) && 
      *(p + offset(nx, ny, x - 1, y + 1)) == 128) 
    watershed(p, nx, ny, x - 1, y + 1, pOut);
  // Bounding raster cell located to the W
  if (inRaster(nx, ny, x - 1, y) && 
      //    !*(pOut + offset(nx, ny, x - 1, y)) && 
      *(p + offset(nx, ny, x - 1, y)) == 1) 
    watershed(p, nx, ny, x - 1, y, pOut);
  // Bounding raster cell located to the NW
  if (inRaster(nx, ny, x - 1, y - 1) && 
      //    !*(pOut + offset(nx, ny, x - 1, y - 1)) && 
      *(p + offset(nx, ny, x - 1, y - 1)) == 2) 
    watershed(p, nx, ny, x - 1, y - 1, pOut);
  // Bounding raster cell located to the N
  if (inRaster(nx, ny, x, y - 1) && 
      //    !*(pOut + offset(nx, ny, x, y - 1)) && 
      *(p + offset(nx, ny, x, y - 1)) == 4) 
    watershed(p, nx, ny, x, y - 1, pOut);
  // Bounding raster cell located to the NE 
  if (inRaster(nx, ny, x + 1, y - 1) && 
      //    !*(pOut + offset(nx, ny, x + 1, y - 1)) && 
      *(p + offset(nx, ny, x + 1, y - 1)) == 8) 
    watershed(p, nx, ny, x + 1, y - 1, pOut);
  
}


//  Function: watershed (version 1, see below)
//
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Cells to be processed are dinamically managed in a queue, up to the basin boundaries
//
//void watershed_v1(int* p, int nx, int ny, int x, int y, int* pOut)
//  void watershed_v1(double* p, int nx, int ny, int pp_offset, int* pOut)
  void watershed_v1(double* p, int nx, int ny, int pp_offset, double* pOut)
{
  int q[10000];    // Queue of raster cells (offset in memory) to be processed
  int delta;      // Offset in memory from base queue address of a raster cell
  int n = 0;      // Number of raster cells to be processed in queue
  int nLoop = 0;  // Counter for loops over cells
  int x,y;
  ///printf("DEBUG: col=%d,row=%d\n",x,y);
  
  // printf("TEST Row, cell 60: %d,%d\n", getRow(nx,ny,60), getCol(nx,ny,60));
  // Set raster cell in the output file
  delta = pp_offset; //offset(nx, ny, x, y);
  *(pOut + delta) = 1;
  *(p + delta) = -10; // EC 20210316
  // Store raster cell offset in the queue  and update number of elements
  q[0] = delta;
  n++;
  
  ///printf("BEFORE n=%d and size(n)=%d\n", n, (int)sizeof(n));
  
  // Process all pending (until any) raster cells in the queue
  //for (int i = 0; i < n; i++) {
  while (n > 0) {
    //printf("DEBUG: IN THE LOOP n=%d\n", n); // REMOVE PRINTF ?? 
    nLoop++;
   // if (nLoop % 10000 == 0) printf("%d ", nLoop);  // Print number of internal loops
    
    // Pick up top raster cell
    x = getCol(nx, ny, q[0]);   // ATTENTION: base 0 or 1?
    y = getRow(nx, ny, q[0]);
    //printf("col=%d row=%d\n", x, y);
    
    // Queue just full. This should be better managed incrementing its size dinamically
    if (n > 9990) {
    //THIS MUST BE MODIFIED  printf("\nAborted! Internal buffer for cells to be processed just full! Size to be incremented!\n");
      return;
    }
    
    // Investigate D8 raster cells all around the cell under investigation
    // Bounding raster cell located to the E
    if (inRaster(nx, ny, x + 1, y) && *(p + offset(nx, ny, x + 1, y)) == 16) {
      delta = offset(nx, ny, x + 1, y);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;                    
    }
    // Bounding raster cell located to the SE
    if (inRaster(nx, ny, x + 1, y + 1) && *(p + offset(nx, ny, x + 1, y + 1)) == 32) {
      delta = offset(nx, ny, x + 1, y + 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;                   
    }
    // Bounding raster cell located to the S
    if (inRaster(nx, ny, x, y + 1) && *(p + offset(nx, ny, x, y + 1)) == 64) {
      delta = offset(nx, ny, x, y + 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the SW
    if (inRaster(nx, ny, x - 1, y + 1) && *(p + offset(nx, ny, x - 1, y + 1)) == 128) {
      delta = offset(nx, ny, x - 1, y + 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the W
    if (inRaster(nx, ny, x - 1, y) && *(p + offset(nx, ny, x - 1, y)) == 1) {
      delta = offset(nx, ny, x - 1, y);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the NW
    if (inRaster(nx, ny, x - 1, y - 1) && *(p + offset(nx, ny, x - 1, y - 1)) == 2) {
      delta = offset(nx, ny, x - 1, y - 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the N
    if (inRaster(nx, ny, x, y - 1) && *(p + offset(nx, ny, x, y - 1)) == 4) {
      delta = offset(nx, ny, x, y - 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the NE 
    if (inRaster(nx, ny, x + 1, y - 1) && *(p + offset(nx, ny, x + 1, y - 1)) == 8) {
      delta = offset(nx, ny, x + 1, y - 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    
    //printf("DEBUG AT THE END - n=%d\n",n);
    
    // Just completed the analysis on current raster cell. 
    // Raster cells in queue are shifted up of one position, current raster cell being
    // the facto removed from the queue, and the number of elements is updated
    for (int i = 0; i < n; i++) q[i] = q[i + 1];
    n--;
  }
}







//
// Function:   resizeQueue
// Scope:      resize (doubling) an existing queue to store integer offsets of raster cells
// Parameters: int* q pointer to existing queue
//             int n  current number of elements
// Return a pointer to the resized queue, while preserving previous values
//
int* resizeQueue(int* q, int n)
{
  int* tmp = (int*)CPLMalloc(2*n*sizeof(int));
  
  printf("resizeQueue function: %d\n", n);
  
  // Copy input queue to the new one element by element. Not initialized elements
  // in the second half of the queue do not need any further action at this stage
  for (int i = 0; i < n; i++) *(tmp + i) = *(q + i);
  // Free up the memory allocated for the previous queue (not needed anymore!);
  CPLFree(q);
  
  return tmp;
}

//
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Current version implements a dinamic queue, up to the basin boundaries
 //
// void watershed_v2(int* p, int nx, int ny, int x, int y, int* pOut)
void watershed_v2(double* p, int nx, int ny, int pp_offset, double* pOut)
{
  int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
  int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
  int delta;        // Offset in memory from base queue address of a raster cell
  int n = 0;        // Number of raster cells to be processed in queue
  int nLoop = 0;    // Counter for loops over cells
  int x,y;
//  printf("DEBUG: col=%d,row=%d\n", x, y);
  
  q = (int*)CPLMalloc(sizeof(int)*qSize);
  
  // printf("TEST Row, cell 60: %d,%d\n", getRow(nx,ny,60), getCol(nx,ny,60));
  // Set raster cell in the output file
  delta = pp_offset; // delta=offset(nx, ny, x, y);
  *(pOut + delta) = 1;
  *(p + delta) = -10; // EC 20210316
  // Store raster cell offset in the queue  and update number of elements
  q[0] = delta;
  n++;
  
  printf("BEFORE n=%d and size(n)=%ld\n", n, sizeof(n));
  
  // Process all pending (until any) raster cells in the queue
  while (n > 0) {
    //printf("DEBUG: IN THE LOOP n=%d\n", n);
    nLoop++;
    if (nLoop % 100000 == 0) printf("%d\n", nLoop);  // Print number of internal loops
    
    //printf("%d\n", n);
    // Pick up top raster cell
    x = getCol(nx, ny, q[0]);   // ATTENTION: base 0 or 1?
    y = getRow(nx, ny, q[0]);
    //printf("col=%d row=%d\n", x, y);
    
    // Queue is just full, only 10 raster cells to accomodate
    if (n > (qSize-10)) {
      q = resizeQueue(q, qSize);
      qSize *= 2;
     // puts("Press any key to continue ... ");
     // getchar();
    }
    
    // Investigate D8 raster cells all around the cell under investigation
    // Bounding raster cell located to the E
    if (inRaster(nx, ny, x + 1, y) && *(p + offset(nx, ny, x + 1, y)) == 16) {
      delta = offset(nx, ny, x + 1, y);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the SE
    if (inRaster(nx, ny, x + 1, y + 1) && *(p + offset(nx, ny, x + 1, y + 1)) == 32) {
      delta = offset(nx, ny, x + 1, y + 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the S
    if (inRaster(nx, ny, x, y + 1) && *(p + offset(nx, ny, x, y + 1)) == 64) {
      delta = offset(nx, ny, x, y + 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the SW
    if (inRaster(nx, ny, x - 1, y + 1) && *(p + offset(nx, ny, x - 1, y + 1)) == 128) {
      delta = offset(nx, ny, x - 1, y + 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the W
    if (inRaster(nx, ny, x - 1, y) && *(p + offset(nx, ny, x - 1, y)) == 1) {
      delta = offset(nx, ny, x - 1, y);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the NW
    if (inRaster(nx, ny, x - 1, y - 1) && *(p + offset(nx, ny, x - 1, y - 1)) == 2) {
      delta = offset(nx, ny, x - 1, y - 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the N
    if (inRaster(nx, ny, x, y - 1) && *(p + offset(nx, ny, x, y - 1)) == 4) {
      delta = offset(nx, ny, x, y - 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    // Bounding raster cell located to the NE
    if (inRaster(nx, ny, x + 1, y - 1) && *(p + offset(nx, ny, x + 1, y - 1)) == 8) {
      delta = offset(nx, ny, x + 1, y - 1);
      *(pOut + delta) = 1;
      q[n] = delta;
      n++;
    }
    
    //printf("DEBUG AT THE END - n=%d\n",n);
    
    // Just completed the analysis on current raster cell.
    // Raster cells in queue are shifted up of one position, current raster cell being
    // the facto removed from the queue, and the number of elements is updated
    for (int i = 0; i < n; i++) q[i] = q[i + 1];
    n--;
  }
  
  CPLFree(q);  // Frees memory allocated to the queue (whether is the original or a resized one)
}























// TO INSERT::: std::vector<double> SpatRaster::readValues(size_t row, size_t nrows, size_t col, size_t ncols){
//Rcpp::IntegerVector SpatRaster::watershed2(int pp_offset,SpatOptions opt) {
SpatRaster  SpatRaster::watershed2(int pp_offset,SpatOptions &opt) {
  // DA TESTARE
  SpatRaster out=geometry();
  //std::vector<std::string> oname="watershed";
  //out.setNames(oname);
  int nx=ncol();
  int ny=nrow();
  //printf("nx=%d ny=%d\n",nx,ny);
  //Rcpp::IntegerVector pOut(nx*ny);
   // https://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027/C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm 
  std::vector<double> p=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
  

 //SEE HERE https://stackoverflow.com/questions/26488480/how-can-i-trust-casting-from-double-to-integer

   //int *q=p.begin();
  //int *qOut=pOut.begin();
  // https://www.google.com/search?q=how+to+express+NumericVector+as+a+pointer&oq=how+to+express+NumericVector+as+a+pointer&aqs=chrome..69i57j33i160.13306j0j15&sourceid=chrome&ie=UTF-8
  // https://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
  // http://adv-r.had.co.nz/Rcpp.html
  
  std::vector<double> pOutv(nx*ny,0);
  // EC 20210319 pOutv.reserve(nx*ny);
  // EC 20210319 std::fill(pOutv.begin(), pOutv.end(), trunc(0));
  

  
  ///see
  //watershed_v1(&p[0],nx,ny,pp_offset,pOut.begin());
  watershed_v2(&p[0],nx,ny,pp_offset,&pOutv[0]);
  if (!out.writeStart(opt)) {
    readStop();
    return out;
  }
  // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
  out.writeValues(pOutv,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
  out.writeStop();
  return out;

  //return(pOut);
  
}

/// 20220809
/// PITFINDER 




// TO INSERT::: std::vector<double> SpatRaster::readValues(size_t row, size_t nrows, size_t col, size_t ncols){
//Rcpp::IntegerVector SpatRaster::watershed2(int pp_offset,SpatOptions opt) {
SpatRaster  SpatRaster::pitfinder2(SpatOptions &opt) {
  // DA TESTARE
  SpatRaster out=geometry();
  //std::vector<std::string> oname="watershed";
  //out.setNames(oname);
  int nx=ncol();
  int ny=nrow();
  //printf("nx=%d ny=%d\n",nx,ny);
  //Rcpp::IntegerVector pOut(nx*ny);
  // https://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027/C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm 
  std::vector<double> p=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
  
  
  //SEE HERE https://stackoverflow.com/questions/26488480/how-can-i-trust-casting-from-double-to-integer
  
  //int *q=p.begin();
  //int *qOut=pOut.begin();
  // https://www.google.com/search?q=how+to+express+NumericVector+as+a+pointer&oq=how+to+express+NumericVector+as+a+pointer&aqs=chrome..69i57j33i160.13306j0j15&sourceid=chrome&ie=UTF-8
  // https://dirk.eddelbuettel.com/code/rcpp/Rcpp-quickref.pdf
  // http://adv-r.had.co.nz/Rcpp.html
  
  std::vector<double> pOutv(nx*ny,0);
  // EC 20210319 pOutv.reserve(nx*ny);
  // EC 20210319 std::fill(pOutv.begin(), pOutv.end(), trunc(0));
  
  
  
  ///see
  pitfinder(&p[0],nx,ny,&pOutv[0]);
  if (!out.writeStart(opt)) {
    readStop();
    return out;
  }
  // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
  out.writeValues(pOutv,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
  out.writeStop();
  return out;
  
  //return(pOut);
  
}

// void watershed_v2(int* p, int nx, int ny, int x, int y, int* pOut)
void pitfinder(double* p, int nx, int ny, double* pOut)
{
  //int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
  //  int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
  //int delta;        // Offset in memory from base queue address of a raster cell
  //int n = 0;        // Number of raster cells to be processed in queue
  //  int nLoop = 0;    // Counter for loops over cells
  int x,y;
  int cnt=1;
  //  int pdown;
  //hu int 
  
  //  printf("DEBUG: col=%d,row=%d\n", x, y);
  
  // ####
  //   
  // ## TROVARE HLOS / BUCHI !!! 
  // #FROM THE CELL(x):  (see help terrein)
  // ##   ## 
  // ## 32	64	128
  // ## 16	x	1
  // ## 8	4	2
  // ##TO THE CELL (y):
  // ## 2 4 8
  // ## 1 y 16
  // ##128   64   32
  //   
  // ####  
  //   
  
  // q = (int*)CPLMalloc(sizeof(int)*qSize);
  
  // printf("TEST Row, cell 60: %d,%d\n", getRow(nx,ny,60), getCol(nx,ny,60));
  // Set raster cell in the output file
  //delta = pp_offset; // delta=offset(nx, ny, x, y);
  //*(pOut + delta) = 1;
  //*(p + delta) = -10; // EC 20210316
  // Store raster cell offset in the queue  and update number of elements
  //q[0] = delta;
  //n++;
  
//  printf("BEFORE n=%d and size(n)=%d\n", n, sizeof(n));
  for (int i = 0; i < nx*ny; i++) {
    *(pOut+i)=0;
    
  }  
  for (int i = 0; i < nx*ny; i++) {
  
   // *(pOut+i)=0; 
   
    // ## 32	64	128
    // ## 16	x	1
    // ## 8	4	2
    
    
    // ## 2 4 8
    // ## 1 y 16
    // ##128   64   32    
    
    x = getCol(nx, ny,i);  
    y = getRow(nx, ny,i);
    printf("\n x=%d ",x);
    printf("y=%d ",y);
    printf("i=%d \n",i);
  //  printf("p=%f ",*(p+i));
  //  printf("cnt=%d ",cnt);
  //  printf("pout=%f ",*(pOut+i));
    if (*(p+i) == 1) {
      // if (inRaster(nx, ny, x + 1, y)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x + 1, y)));
      // }
      if (inRaster(nx, ny, x + 1, y) && *(p+offset(nx, ny, x + 1, y)) == 16) {
        
       
        *(pOut+i)=*(pOut+offset(nx, ny, x + 1, y));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
      // *(pOut+i)=1;   // TO COMMENT          
      } 
    } else if (*(p+i) == 2) {
      // if (inRaster(nx, ny, x + 1, y+1)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x + 1, y+1)));
      // }
      if (inRaster(nx, ny, x + 1, y+1) && *(p+offset(nx, ny, x + 1, y+1)) == 32) {
      
        *(pOut+i)=*(pOut+offset(nx, ny, x + 1, y+1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
      // *(pOut+i)=1;   // TO COMMENT 
      } 
            
    } else if (*(p+i) == 4) {
      // if (inRaster(nx, ny, x , y-1)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x , y+1)));
      // }
      if (inRaster(nx, ny, x, y+1) && *(p+offset(nx, ny, x, y+1)) == 64) {
        //printf("pit p=%f \n",*(p+i));
        *(pOut+i)=*(pOut+offset(nx, ny, x, y+1));
        
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
        cnt++;
        }
        // *(pOut+i)=1;   // TO COMMENT 
      } 
            
    } else if (*(p+i) == 8) {
      // if (inRaster(nx, ny, x - 1, y+1)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x - 1, y+1)));
      // }
      if (inRaster(nx, ny, x-1, y+1) && *(p+offset(nx, ny, x - 1, y+1)) == 128) {
        
        *(pOut+i)=*(pOut+offset(nx, ny, x-1, y+1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
        // *(pOut+i)=1;   // TO COMMENT 
      } 
    } else if (*(p+i) == 16) {
      // if (inRaster(nx, ny, x - 1, y)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x - 1, y)));
      // }
      if (inRaster(nx, ny, x-1, y) && *(p+offset(nx, ny, x - 1, y)) == 1) {
       
        *(pOut+i)=*(pOut+offset(nx, ny, x-1,y));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
        // *(pOut+i)=1;   // TO COMMENT 
      } 
    } else if (*(p+i) == 32) {
      // if (inRaster(nx, ny, x - 1, y-1)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x - 1, y-1)));
      // }
      if (inRaster(nx, ny, x - 1, y-1) && *(p+offset(nx, ny, x-1, y-1)) == 2) {
      
        *(pOut+i)=*(pOut+offset(nx, ny, x-1,y-1));
        if (*(pOut+i)==0) {
         *(pOut+i)=(double)cnt;
         cnt++;
        }
        // *(pOut+i)=1;   // TO COMMENT 
        
      } 
    } else if (*(p+i) == 64) {
      // if (inRaster(nx, ny, x , y-1)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x, y-1)));
      // }
      if (inRaster(nx, ny, x, y-1) && *(p+offset(nx, ny, x, y-1)) == 4) {
       
        *(pOut+i)=*(pOut+offset(nx, ny, x,y-1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
      } 
    } else if (*(p+i) == 128) {
      // if (inRaster(nx, ny, x + 1, y-1)) {
      //   printf("pit p=%f ,,",*(p+i));
      //   printf("pit p=%f ,,",*(p+offset(nx, ny, x + 1, y-1)));
      // }
      if (inRaster(nx, ny, x + 1, y-1) && *(p+offset(nx, ny, x + 1, y-1)) == 8) {
      
        *(pOut+i)=*(pOut+offset(nx, ny, x+1,y+1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
        // *(pOut+i)=1;   // TO COMMENT 
      } 
    }
    
   // printf("pout=%f",*(pOut+i));
   //  printf("cnt=%d \n",cnt);
  }

}


