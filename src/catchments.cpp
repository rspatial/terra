#include "spatRaster.h"

// C/C++ code
// Author: Ezio Crestaz,Emanuele Cordano
// Date: February 2021 - July 2026
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer 1D array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function recursively call itself up to the watershed boundary or raster limit

// Function: offset
// Scope: return offset of a raster cell stored in a 1D integers array with respect to base address
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
int offset(int nx, int ny, int x, int y) {
  return y * nx + x; //according to original Ezio's code BY ROWS 
}

// Functions: getRow, getCol
// Scope: return offset of a raster cell stored in a 1D integers array with respect to base address
// nx,ny: number of cells in x (longitude) and in y (latitude)
// offset: 
// NOTE: it assumes the offset is a valid value
int getRow(int nx, int ny, int offset) {
  return offset/nx;// according to Ezio's original code // BY ROWS
}

int getCol(int nx, int ny, int offset) {
// according to Ezio's original code
  return offset % nx;// according to Ezio's original code // BY ROWS
}


//  Function: inRaster
// Scope: check if a given cell is valid or not, that's within the raster
//  nx,ny: number of cells in x (longitude) and in y (latitude)
//  x,y: indexes of cell upstream of which the watershed must be computed
bool inRaster(int nx, int ny, int x, int y) {
  return !(x < 0 || x >= nx || y >= ny || y < 0);
}


//  Function: watershed (version 0, recursive)
void watershed(double* p, int nx, int ny, int x, int y, int* pOut) {
//  static int nCall = 0;
//  nCall++;
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
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Cells to be processed are dinamically managed in a queue, up to the basin boundaries
//
//void watershed_v1(int* p, int nx, int ny, int x, int y, int* pOut)
//  void watershed_v1(double* p, int nx, int ny, int pp_offset, int* pOut)
void watershed_v1(double* p, int nx, int ny, int pp_offset, double* pOut) {
  int q[10000];    // Queue of raster cells (offset in memory) to be processed
  int delta;      // Offset in memory from base queue address of a raster cell
  int n = 0;      // Number of raster cells to be processed in queue
//  int nLoop = 0;  // Counter for loops over cells
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
   // nLoop++;
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



// Function:   resizeQueue
// Scope:      resize (doubling) an existing queue to store integer offsets of raster cells
// Parameters: int* q pointer to existing queue
//             int n  current number of elements
// Return a pointer to the resized queue, while preserving previous values
//
int* resizeQueue(int* q, int n)
{
  int* tmp = (int*)CPLMalloc(2*n*sizeof(int));

  //printf("resizeQueue function: %d\n", n); //EC 20231129

  // Copy input queue to the new one element by element. Not initialized elements
  // in the second half of the queue do not need any further action at this stage
  for (int i = 0; i < n; i++) *(tmp + i) = *(q + i);
  // Free up the memory allocated for the previous queue (not needed anymore!);
  CPLFree(q);

  return tmp;
}


// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Current version implements a dinamic queue, up to the basin boundaries
 //
// void watershed_v2(int* p, int nx, int ny, int x, int y, int* pOut)
void watershed_v2(double* p, int nx, int ny, int pp_offset, double* pOut) {
  int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
  int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
  int delta;        // Offset in memory from base queue address of a raster cell
  int n = 0;        // Number of raster cells to be processed in queue
//  int nLoop = 0;    // Counter for loops over cells
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

 // commented by EC 20231129 printf("BEFORE n=%d and size(n)=%ld\n", n, sizeof(n));

  // Process all pending (until any) raster cells in the queue
  while (n > 0) {
    //printf("DEBUG: IN THE LOOP n=%d\n", n);
 //   nLoop++;
   // EC 20231129 if (nLoop % 100000 == 0) printf("%d\n", nLoop);  // Print number of internal loops

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
  if (!out.writeStart(opt,filenames())) {
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

// void watershed_v2(int* p, int nx, int ny, int x, int y, int* pOut)
void pitfinder(double* p, int nx, int ny, double* pOut,int pits_on_boundary)  {

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

//  printf("BEFORE n=%d and size(n)=%d\n", n, sizeof(n));
  for (int i = 0; i < nx*ny; i++) {

    if (pits_on_boundary==0){
      *(pOut+i)=*(p+i)*0; //20241106
    } else {
      *(pOut+i)=0;
    }


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
    if (*(p+i) == 1) {
      if (inRaster(nx, ny, x + 1, y) && *(p+offset(nx, ny, x + 1, y)) == 16) {
        *(pOut+i)=*(pOut+offset(nx, ny, x + 1, y));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;       
        }
      } 
    } else if (*(p+i) == 2) {
      if (inRaster(nx, ny, x + 1, y+1) && *(p+offset(nx, ny, x + 1, y+1)) == 32) {      
        *(pOut+i)=*(pOut+offset(nx, ny, x + 1, y+1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
      }             
    } else if (*(p+i) == 4) {
      if (inRaster(nx, ny, x, y+1) && *(p+offset(nx, ny, x, y+1)) == 64) {
        *(pOut+i)=*(pOut+offset(nx, ny, x, y+1));       
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
       }
      } 

    } else if (*(p+i) == 8) {
      if (inRaster(nx, ny, x-1, y+1) && *(p+offset(nx, ny, x - 1, y+1)) == 128) {

        *(pOut+i)=*(pOut+offset(nx, ny, x-1, y+1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;      
        }
      } 
    } else if (*(p+i) == 16) {
      if (inRaster(nx, ny, x-1, y) && *(p+offset(nx, ny, x - 1, y)) == 1) {       
        *(pOut+i)=*(pOut+offset(nx, ny, x-1,y));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
      } 
    } else if (*(p+i) == 32) {
      if (inRaster(nx, ny, x - 1, y-1) && *(p+offset(nx, ny, x-1, y-1)) == 2) {      
        *(pOut+i)=*(pOut+offset(nx, ny, x-1,y-1));
        if (*(pOut+i)==0) {
         *(pOut+i)=(double)cnt;
         cnt++;
        }

      } 
    } else if (*(p+i) == 64) {
      if (inRaster(nx, ny, x, y-1) && *(p+offset(nx, ny, x, y-1)) == 4) {

        *(pOut+i)=*(pOut+offset(nx, ny, x,y-1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;

        } 

      } 
    } else if (*(p+i) == 128) {
      if (inRaster(nx, ny, x + 1, y-1) && *(p+offset(nx, ny, x + 1, y-1)) == 8) {      
        *(pOut+i)=*(pOut+offset(nx, ny, x+1,y+1));
        if (*(pOut+i)==0) {
          *(pOut+i)=(double)cnt;
          cnt++;
        }
      } 
    } else if (*(p+i)==0){

         if (*(pOut+i)==0 && inRaster(nx, ny, x + 1, y+1)) *(pOut+i)=*(pOut+offset(nx, ny, x+1,y+1));
         if (*(pOut+i)==0 && inRaster(nx, ny, x , y+1)) *(pOut+i)=*(pOut+offset(nx, ny, x,y+1));
         if (*(pOut+i)==0 && inRaster(nx, ny, x - 1, y+1)) *(pOut+i)=*(pOut+offset(nx, ny, x-1,y+1));

         if (*(pOut+i)==0 && inRaster(nx, ny, x + 1, y)) *(pOut+i)=*(pOut+offset(nx, ny, x+1,y));
         if (*(pOut+i)==0 && inRaster(nx, ny, x , y)) *(pOut+i)=*(pOut+offset(nx, ny, x,y));
         if (*(pOut+i)==0 && inRaster(nx, ny, x - 1, y)) *(pOut+i)=*(pOut+offset(nx, ny, x-1,y));

         if (*(pOut+i)==0 && inRaster(nx, ny, x + 1, y-1)) *(pOut+i)=*(pOut+offset(nx, ny, x+1,y-1));
         if (*(pOut+i)==0 && inRaster(nx, ny, x , y-1)) *(pOut+i)=*(pOut+offset(nx, ny, x,y-1));
         if (*(pOut+i)==0 && inRaster(nx, ny, x - 1, y-1)) *(pOut+i)=*(pOut+offset(nx, ny, x-1,y-1));
         if (*(pOut+i)==0) {
           *(pOut+i)=(double)cnt;
           cnt++;
         } 
      }
  }
}



// TO INSERT::: std::vector<double> SpatRaster::readValues(size_t row, size_t nrows, size_t col, size_t ncols){
//Rcpp::IntegerVector SpatRaster::watershed2(int pp_offset,SpatOptions opt) {
SpatRaster  SpatRaster::pitfinder2(int pits_on_boundary, SpatOptions &opt) {
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
  pitfinder(&p[0],nx,ny,&pOutv[0],pits_on_boundary);
  if (!out.writeStart(opt,filenames())) {
    readStop();
    return out;
  }
  // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
  out.writeValues(pOutv,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
  out.writeStop();
  return out;

  //return(pOut);

}

///////////////////////////////////////
//// flow accumulation functions //////
///////////////////////////////////////
void NextCell(double* p, int nx, int ny,int* pnext) {

    //int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
    //  int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
    //int delta;        // Offset in memory from base queue address of a raster cell
    //int n = 0;        // Number of raster cells to be processed in queue
    //  int nLoop = 0;    // Counter for loops over cells
    int i;

    // ## 32	64	128
    // ## 16	x	1
    // ## 8	4	2

    for (i=0;i<nx*ny;i++) {
      *(pnext+i)=i; //-9999;
    }
    for (int x=0;x<nx;x++) {
      for (int y=0;y<ny;y++) {

        i = offset(nx,ny,x,y);
        if (inRaster(nx, ny, x + 1, y) && *(p+i)==1) {
           *(pnext+i)=offset(nx, ny, x + 1, y);
        } else if (inRaster(nx, ny, x+1,y+1) && *(p+i)==2) {

          *(pnext+i)=offset(nx, ny, x+1, y+1);
        } else if (inRaster(nx, ny, x,y+1) && *(p+i)==4) {
          *(pnext+i)=offset(nx, ny, x, y+1);
        } else if (inRaster(nx, ny, x-1,y+1) && *(p+i)==8) {
          *(pnext+i)=offset(nx, ny, x-1, y+1);
        } else if (inRaster(nx, ny, x-1,y) && *(p+i)==16) {
          *(pnext+i)=offset(nx, ny, x-1, y);
        } else if (inRaster(nx, ny, x-1,y-1) && *(p+i)==32) {
          *(pnext+i)=offset(nx, ny, x-1, y-1);
        } else if (inRaster(nx, ny, x,y-1) && *(p+i)==64) {
          *(pnext+i)=offset(nx, ny, x, y-1);
        } else if (inRaster(nx, ny, x+1,y-1) && *(p+i)==128) {
          *(pnext+i)=offset(nx, ny, x+1, y-1);
      } 
    }

  }  

}  


// NIPD 
void NIDP(int* pnext, int nx, int ny,double* nidp_value) {

  int i=0,pp=0;
  int x,y;
  double cnt0; 

  for (i=0;i<nx*ny;i++) {

    *(nidp_value+i)=0;
  }  

  for(x=0;x<nx;x++) {
    for(y=0;y<ny;y++) {

      i = offset(nx,ny,x,y);
      pp=*(pnext+i);
      if (pp!=-9999) {
        cnt0=*(nidp_value+pp);
        cnt0++;
        *(nidp_value+pp)=cnt0;
        cnt0=0;
      }    
    }
  }  
}  


 // TO INSERT::: std::vector<double> SpatRaster::readValues(size_t row, size_t nrows, size_t col, size_t ncols){
 //Rcpp::IntegerVector SpatRaster::watershed2(int pp_offset,SpatOptions opt) {
 SpatRaster  SpatRaster::NIDP2(SpatOptions &opt) {
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


   std::vector<double> pOutv(nx*ny,0);
   std::vector<int> pnext(nx*ny,0);
   std::vector<double> nidp_value(nx*ny,0);

   NextCell(&p[0],nx,ny,&pnext[0]);
   NIDP(&pnext[0],nx,ny,&nidp_value[0]); 
   if (!out.writeStart(opt,filenames())) {
     readStop();
     return out;
   }
   // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
   out.writeValues(nidp_value,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
   out.writeStop();
   return out;

   //return(pOut);

 }

// FlowAccu algorithm 5 
// Reference: https://link.springer.com/article/10.1007/s11707-018-0725-9
void FlowAccu(int* pnext, int nx, int ny,double* nidp_value,double* flowaccu_value) {

  int n=0,flowpath_cond;
//  int x,y;
  double nAccu=0; 

  for (int i=0;i<nx*ny;i++) {

    *(flowaccu_value+i)=1;

  }  
  for (int i=0;i<nx*ny;i++) if (*(nidp_value+i)==0) {
    //printf("i=%d \n",i);
    //printf("flowpath_cond=%d \n",flowpath_cond);
    n=i;
    nAccu=0;
    flowpath_cond=1;
    //printf("n=%d \n",n);
    //printf("flowpath_cond=%d \n",flowpath_cond);
    do {
      *(flowaccu_value+n)+=nAccu;
      nAccu=*(flowaccu_value+n);
      if (*(nidp_value+n)>=2) {
        *(nidp_value+n)-=1;
        flowpath_cond=0;
      } else {
        n=*(pnext+n);
      }
      //printf("n=%d \n",n);
      //printf("flowpath_cond=%d \n",flowpath_cond);
    } while (flowpath_cond==1);
  }
}  


// It is added an array of weight to each cell (e.g. area, averaged precipitation)
void FlowAccu_weight(int* pnext, int nx, int ny,double* nidp_value,double* flowaccu_value,double* weight) {

  int n=0, flowpath_cond;
  double nAccu=0; 

  for (int i=0;i<nx*ny;i++) {

    *(flowaccu_value+i)=*(weight+i);

  }  
  for (int i=0;i<nx*ny;i++) if (*(nidp_value+i)==0) {
    //printf("i=%d \n",i);
    //printf("flowpath_cond=%d \n",flowpath_cond);
    n=i;
    nAccu=0;
    flowpath_cond=1;
    //printf("n=%d \n",n);
    //printf("flowpath_cond=%d \n",flowpath_cond);
    do {
      *(flowaccu_value+n)+=nAccu;
      nAccu=*(flowaccu_value+n);
      if (*(nidp_value+n)>=2) {
        *(nidp_value+n)-=1;
        flowpath_cond=0;
      } else {
        n=*(pnext+n);
      }
      //printf("n=%d \n",n);
      //printf("flowpath_cond=%d \n",flowpath_cond);
    } while (flowpath_cond==1);
  }
}  



SpatRaster  SpatRaster::flowAccu2(SpatOptions &opt) {
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


  std::vector<double> pOutv(nx*ny,0);
  std::vector<int> pnext(nx*ny,0);
  std::vector<double> nidp_value(nx*ny,0);
  std::vector<double> flowaccu_value(nx*ny,1);


  NextCell(&p[0],nx,ny,&pnext[0]);
  NIDP(&pnext[0],nx,ny,&nidp_value[0]); 
  FlowAccu(&pnext[0],nx,ny,&nidp_value[0],&flowaccu_value[0]);
  if (!out.writeStart(opt,filenames())) {
    readStop();
    return out;
  }
  // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
  out.writeValues(flowaccu_value,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
  out.writeStop();
  return out;

  //return(pOut);

}

SpatRaster  SpatRaster::flowAccu2_weight(SpatRaster weight,SpatOptions &opt) {
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
  std::vector<double> weigh=weight.getValues(0,opt);

  std::vector<double> pOutv(nx*ny,0);
  std::vector<int> pnext(nx*ny,0);
  std::vector<double> nidp_value(nx*ny,0);
  std::vector<double> flowaccu_value(nx*ny,1);


  NextCell(&p[0],nx,ny,&pnext[0]);
  NIDP(&pnext[0],nx,ny,&nidp_value[0]); 
  FlowAccu_weight(&pnext[0],nx,ny,&nidp_value[0],&flowaccu_value[0],&weigh[0]);

  if (!out.writeStart(opt,filenames())) {
    readStop();
    return out;
  }
  // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
  out.writeValues(flowaccu_value,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
  out.writeStop();
  return out;


}

/////// FLOW DIRECTION 

// #define NODATA_2 -999999

// ESRI terra covension:
//   //     # x-1 x 	x+1
//   // y-1 ## 32	64	128
//   // y   ## 16	0	1
//   // y+1 ## 8	4	2
//   

//  Li at al, 2022's convection 
//   //     # x-1 x 	x+1
//   // y-1 ## 1	2	3
//   // y   ## 8	0	4
//   // y+1 ## 7	6	5

/* nx,ny number of cols,rows 
 * x,y    row,col
 * pv flow direction value 
 * conv_type convention type used , see code */

int nextcell_point_conv1(int nx, int ny,int x,int y,double pv,int conv_type) {

    std::vector<double> idir={0,1,2,4,8,16,32,64,128}; // ESRI conv: 0 cell, ... (clockwise from right) 
    int iout=offset(nx,ny,x, y);
    if (conv_type==2) {
     // Li at al 2021 (from Orlandini at al 20)
     idir={0,4,5,6,7,8,1,2,3}; // Li-Orlandini conv 0 cell, ... (clockwise from right) 
    }

    if (inRaster(nx, ny, x + 1, y) && pv==idir[1]) {
        iout=offset(nx, ny, x + 1, y);
    } else if (inRaster(nx, ny, x+1,y+1) && pv==idir[2]) {

      iout=offset(nx, ny, x+1, y+1);
    } else if (inRaster(nx, ny, x,y+1) && pv==idir[3]) {
      iout=offset(nx, ny, x, y+1);
    } else if (inRaster(nx, ny, x-1,y+1) && pv==idir[4]) {
      iout=offset(nx, ny, x-1, y+1);
    } else if (inRaster(nx, ny, x-1,y) && pv==idir[5]) {
      iout=offset(nx, ny, x-1, y);
    } else if (inRaster(nx, ny, x-1,y-1) && pv==idir[6]) {
      iout=offset(nx, ny, x-1, y-1);
    } else if (inRaster(nx, ny, x,y-1) && pv==idir[7]) {
      iout=offset(nx, ny, x, y-1);
    } else if (inRaster(nx, ny, x+1,y-1) && pv==idir[8]) {
      iout=offset(nx, ny, x+1, y-1);
    } else {
      iout=offset(nx, ny, x, y);
    }
    //Rprintf("iout=%d\n",iout);
    return(iout);
    }

  /* e elevation
   * nx,ny number of cols,rows
   * sr SR in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * sm SM in Li et al,2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * sfacet drainage facet 
   * tdc TD_c in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * tdd TD_c in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   */

void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,
                     double* tdc, double *tdd,double L,
                     std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,int conv_type,int use_lad) {

  //int i;
  ///int L=1; // dx=dy=L=1
  //int ncell=nx*ny;
  double x,y;
 // int efacet=0;
  double e0,e1,e2;
  //int facet=0;
//  std::vector<double> ddp1 = {2,2,4,4,6,6,8,8}; // this routine uses Orlandini-Li et, 2022's  convention
//  std::vector<double> ddp2 = {1,3,3,5,5,7,7,1};
// std::vector<double> sigma = {0,1,-1,1,-1,1,-1,1,-1};
 // int nncell=8;
  //std::vector<double> vOut(nx*ny,0);
  //std::vector<double> slope_mgn(ncell,0);
 // double slope=
   for (int i = 0; i < nx*ny; i++) {

    e0=*(e+i);
    std::vector<int> offset1(nncell,i);
    std::vector<int> offset2(nncell,i);
    x = getCol(nx, ny, i);   // ATTENTION: base 0 or 1?
    y = getRow(nx, ny, i);

    double slope_mgn=0; //-1;
    double slope_mgn_temp=slope_mgn;
    double mean_e=e0;// not 0 corrected on 20260109
    double mean_e_temp=mean_e;// not 0 corrected on 20260109
    int facet=0;
    //double flow_angle_tan=0;
    for (int j=0; j<nncell; j++){

      //if (ddp1[j]==1) {

    e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[j],conv_type));
    e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[j],conv_type));
    //Rprintf("ba2 j=%d facet=%d i=%d e0=%f e1=%f e2=%f\n",j,facet,i,e0,e1,e2);
    //if (std::isnan(e1) | std::isnan(e0) | e1==e0) e1=1.05*e0; 
    //if (std::isnan(e2) | std::isnan(e0) | e2==e0) e2=1.05*e0;

    /// facet=j;
    ///double mean_e_temp=(e0+e1+e2)/3; //EXPERIMENTAL    /////pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L;
    mean_e_temp=(e0+e1+e2)/3;
    if ((e0==e1) & (e1==e2)) {

      slope_mgn_temp=-1; // added on 20240429
    } else if ((e0>=e1) & (e1>=e2)) { // 210
      double flow_angle_tan_temp=0;
      if (e0!=e1) flow_angle_tan_temp=(e1-e2)/(e0-e1); //ec 20251023
      if (flow_angle_tan_temp<=1) {
        slope_mgn_temp=pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L; //210
      } else {

        slope_mgn_temp=((e0-e2)/(L*sqrt(2)));//210
      }

    //  slope_mgn_temp1=((e0-e1)/L); // 120 102     
     // double flow_angle_tan=((e1-e2)/(e0-e1));
    } else if ((e0>=e1) & (e1<=e2)) { // 120 102
      slope_mgn_temp=((e0-e1)/L); // 120 102
     // double flow_angle_tan=0;
    } else if ((e0>=e2)  & (e1>=e0)) { //201
      slope_mgn_temp=((e0-e2)/(L*sqrt(2)));
    }  else { 
      slope_mgn_temp=-1;
      // 012 021
      //nothing 
    }

    if (slope_mgn_temp>slope_mgn) {
        slope_mgn=slope_mgn_temp;
        mean_e=mean_e_temp;
        facet=j;        
    } else  if (slope_mgn_temp==slope_mgn) {
         if (mean_e_temp<mean_e) {
          slope_mgn=slope_mgn_temp;
          mean_e=mean_e_temp;
          facet=j;
          } else if (mean_e_temp==mean_e){
              facet=0; // uncommented on 20250109 // more actions to do // 2026 01 18           
          }
     }
    }

    double e_temp=e0;
    if (facet==0)  for (int j=0; j<nncell; j++){
      e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[j],conv_type));
      e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[j],conv_type));
      mean_e_temp=(e0+e1+e2)/3;
      if (!(std::isnan(e1) || std::isnan(e2))) { // to rewrite 
        if (e1<e_temp) {
          e_temp=e1;
          mean_e=mean_e_temp;
          facet=j;
        } 
        if (e2<e_temp) {
          e_temp=e2;
          mean_e=mean_e_temp;
          facet=j;          
        } // added on 2024 09 16
      }
    }
    e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[facet],conv_type));
    e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[facet],conv_type)); 
    //if (std::isnan(e1) | std::isnan(e0) | e1==e0) e1=1.05*e0;
    //if (std::isnan(e2) | std::isnan(e0) | e2==e0) e2=1.05*e0;
    *(sfacet+i)=facet;
    //*(sr+i)=(std::atan((e1-e2)/(e0-e1))); // RAD 20240912
    if (e1==e0) {
      if (e1<=e2) {
        *(sr+i)=0;
      } else {
        *(sr+i)=M_PI/4;
      }
    } else {
      *(sr+i)=(std::atan((e1-e2)/(e0-e1))); // mod ec 20250322 //RAD 20240912
    }
    *(sm+i)=pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L;
    if (*(sr+i)<0) {

      *(sr+i)=0;
      *(sm+i)=abs(e0-e1)/L;

    } else if (*(sr+i)>M_PI/4) {

        *(sr+i)=M_PI/4;
        *(sm+i)=abs(e0-e2)/(L*sqrt(2));
    }

    if (use_lad==1) {
      *(tdc+i)=(*(sr+i));
      *(tdd+i)=(M_PI/4-*(sr+i));
    } else {
      *(tdc+i)=L*std::sin(*(sr+i));
      *(tdd+i)=sqrt(2)*L*std::sin(M_PI/4-*(sr+i));
    }
  }
}


// returns false if the maximum number of iterations was exceeded
bool transverse_deviation(double *e, double *tdc, double *tdd,double *sr,double *sm, int *sfacet,int nx, int ny, double L,
		double *atdc, double *atdd, double *atdplus,double *atdplus0, double *pflow,int *has_upstream,int *kupdate,
		double *nidps, double lambda,  std::vector<double> ddp1, std::vector<double> ddp2, std::vector<double> sigma,
		int nncell, int conv_type, int use_lad, int max_iters) {   

  bool ok = true;
  int x,y;

  //int k=1;
  //int niter=nx*ny;
  int exit_cond=0;
  int exit_cond1=0;
  int facet; //,nextc,nextd;
  int nextp=0;
  double atdplus_temp=0;
  double e0,e1,e2;
  double pflow_estimate=ddp1[0];

  for (int i = 0; i < nx*ny; i++) {
    *(has_upstream+i)=0;
  } 

  for (int i = 0; i < nx*ny; i++) { 
   x = getCol(nx, ny, i);   // ATTENTION: base 0 or 1?
   y = getRow(nx, ny, i); 
   facet=*(sfacet+i);
   ////
   e0=*(e+i);
   int nextc=nextcell_point_conv1(nx,ny,x,y,ddp1[facet],conv_type);
   int nextd=nextcell_point_conv1(nx,ny,x,y,ddp2[facet],conv_type);
   e1=*(e+nextc);
   e2=*(e+nextd); 

   ////



   *(atdc+i)=*(tdc+i)*sigma[facet];
   *(atdplus+i)=0;
   *(atdd+i)=*(tdd+i)*sigma[facet]*(-1); // 20240912
   *(pflow+i)=ddp1[0]; // initialization of pflow estimate .

   *(kupdate+i)=0;
   if ((e0<=e1) & (e0>e2)) { //// pit ???  //// coorect here 20250116
     *(pflow+i)=ddp2[facet]; // initialization of pflow estimate .
     *(atdplus+i)= *(atdd+i);
     //*(pflow2+i)=ddp1[facet]; //tocontinue...
   } else if ((e0<=e2) & (e0>e1)) {
     *(pflow+i)=ddp1[facet]; // initialization of pflow estimate .
     *(atdplus+i)= *(atdc+i);
   } else if ((abs(*(atdc+i))<abs(*(atdd+i))) || ((abs(*(atdc+i))<=abs(*(atdd+i))) & (use_lad==1))){
     *(pflow+i)=ddp1[facet]; // initialization of pflow estimate .
     *(atdplus+i)= *(atdc+i);  
   } else if (abs(*(atdc+i))>=abs(*(atdd+i))) {
     *(pflow+i)=ddp2[facet]; // initialization of pflow estimate .
     *(atdplus+i)= *(atdd+i);

   }
   *(atdplus0+i)=*(atdplus+i);
   pflow_estimate=*(pflow+i);
   nextp=nextcell_point_conv1(nx,ny,x,y,pflow_estimate,conv_type);
   if (nextp!=i) {

     *(has_upstream+nextp)=1+*(has_upstream+nextp);
   }
   *(kupdate+i)=0;
  } // intialization of atdplus+i)

  int cnt=0;
  int cnt1=0;
 if (lambda>0) do { 
  for (int i = 0; i < nx*ny; i++) {
     *(has_upstream+i)=0;
  }
  for (int i = 0; i < nx*ny; i++) {
     x = getCol(nx, ny, i);   // ATTENTION: base 0 or 1?
     y = getRow(nx, ny, i); 
     pflow_estimate=*(pflow+i);
     nextp=nextcell_point_conv1(nx,ny,x,y,pflow_estimate,conv_type);
     if (nextp!=i) {
        *(has_upstream+nextp)=1+*(has_upstream+nextp);
     }
  }

  for (int j = 0; j < nx*ny; j++) {
  int i=j; 

  if ((*(kupdate+j)==0) & ((*(has_upstream+j)==0) & (cnt1>=0))) cnt++; // ???
  if ((*(kupdate+j)==0) & ((*(has_upstream+j)==0) & (cnt1>=0))) do {
    *(atdplus+j)=*(atdplus0+j); // ?????
   /// *(has_upstream+j)=-2;
    exit_cond=0;
    *(kupdate+i)=cnt;
    x = getCol(nx, ny, i);   // ATTENTION: base 0 or 1
    y = getRow(nx, ny, i);
   // nextp=i;
  //  pflow_estimate=*(pflow+j);
   // pflow_estimate0=*(pflow+i);
   // e0=*(e+i); // not j ec 20260430

    nextp=nextcell_point_conv1(nx,ny,x,y,*(pflow+i),conv_type);

    // analyse nextp cell 
    int xp = getCol(nx, ny, nextp);   // ATTENTION: base 0 or 1
    int yp = getRow(nx, ny, nextp);

    facet=*(sfacet+nextp); //*(sfacet+i);
    e0=*(e+nextp);

  //  int pflow_estimate0=*(pflow+nextp);
  //  int nextp0=nextcell_point_conv1(nx,ny,xp,yp,pflow_estimate0,conv_type);
    int nextpc=nextcell_point_conv1(nx,ny,xp,yp,ddp1[facet],conv_type);
    int nextpd=nextcell_point_conv1(nx,ny,xp,yp,ddp2[facet],conv_type);
    e1=*(e+nextpc);
    e2=*(e+nextpd); 
    pflow_estimate=*(pflow+nextp);
    int nextq=nextcell_point_conv1(nx,ny,x,y,pflow_estimate,conv_type); // 20250430
    double atdplus_nextpc=*(tdc+nextp)*sigma[facet]+*(atdplus+i)*lambda;
    double atdplus_nextpd=*(tdd+nextp)*sigma[facet]*(-1)+*(atdplus+i)*lambda; // 20240924
    if ((e0<=e1) & (e0>e2)) { //// pit ???  //// coorect here 20250116
        //*(atdplus+nextd)=*(atdplus+nextd)+*(atdd+i); //DA VEDERE BENE
      pflow_estimate=ddp2[facet];
      nextq=nextpd;
        //*(atdplus+nextp)
      atdplus_temp=atdplus_nextpd;

        //*(pflow2+i)=ddp1[facet]; //tocontinue...
    } else if ((e0<=e2) & (e0>e1)) {

      nextq=nextpc;// ERRORE!!!
      pflow_estimate=ddp1[facet];
      atdplus_temp=atdplus_nextpc;
    } else if ((abs(atdplus_nextpc)<abs(atdplus_nextpd)) || ((abs(atdplus_nextpc)<=abs(atdplus_nextpd)) & (use_lad==1))){

      pflow_estimate=ddp1[facet];
      nextq=nextpc; // cardinal
      atdplus_temp=atdplus_nextpc;
    } else if (abs(atdplus_nextpd)<=abs(atdplus_nextpc)) {
   // } else  {
      pflow_estimate=ddp2[facet];
      nextq=nextpd; // diagonal
      atdplus_temp=atdplus_nextpd;

    } else {

      double slo1=(e0-e1)/L;
      double slo2=(e0-e2)/(L*sqrt(2));
     // Rprintf("HERE!!! x=%d y=%d i=%d slo1=%f slo2=%f e0=%f e1=%f e2=%f L=%f\n",x,y,i,slo1,slo2,e0,e1,e2,L);
      if (slo2>=slo1) {
        pflow_estimate=ddp2[facet];
        nextq=nextpd;
        atdplus_temp=atdplus_nextpd;

      } else if (slo2<slo1) {
        pflow_estimate=ddp1[facet];
        nextq=nextpc;
        atdplus_temp=atdplus_nextpc;
      }
    }
    // work here 20260122
    // condizione sopre ok se non cambia direzione ??? se e uguale e cambia dierezione???
    if ((abs(atdplus_temp)>=abs(*(atdplus+nextp)))) { // 20260119    if (abs(atdplus_temp)>=abs(*(atdplus+i))){ // 20260119
        // ADD A CONTROL (kupdate+nextp )
      *(atdplus+nextp)=atdplus_temp; // corrected on 20260429
      *(pflow+nextp)=pflow_estimate;
      *(kupdate+nextp)=cnt;
    } else {
    }
    if (nextp!=i) { // 20260429      
      i=nextp;
      nextp=nextq;
      exit_cond=0;
    } else {
      exit_cond=2;
    }
  } while (exit_cond==0);
  cnt1++;
  exit_cond1=2;

  // for statemant to verify kupdate!=0
  }
  for (int j = 0; j < nx*ny; j++) {    
    if (*(kupdate+j)==0) exit_cond1=0;
  }
  if (cnt1>max_iters) {
    ok = false;
    exit_cond1=2;       

  }  
 } while  (exit_cond1==0);
  // NOVALUE
  for (int i=0;i<nx*ny;i++) {
    double e0=*(e+i);
    if (std::isnan(e0)) {
      *(pflow+i)=e0;
   }
  }
  return ok;
}


bool d8ltd_computation(double *e,int nx,int ny,double L,double lambda,int use_lad,int max_iters,double*pflow) {
  //std::vector<double> pOutv(nx*ny,0);
  std::vector<int> has_upstream(nx*ny,0);
  std::vector<int> sfacet(nx*ny,0);
  std::vector<double> tdc(nx*ny,0);
  std::vector<double> tdd(nx*ny,0);
  std::vector<double> atdc(nx*ny,0);
  std::vector<double> atdd(nx*ny,0);
  std::vector<double> atdplus(nx*ny,0);
  std::vector<double> atdplus0(nx*ny,0);
  std::vector<double> sr(nx*ny,0);
  std::vector<double> sm(nx*ny,0);


  std::vector<int> kupdate(nx*ny,0);
  std::vector<double> npids(nx*ny,0);  // npid number 
  /* FLOW Direction Convention */ 

  //std::vector<double> ddp1 = {0,2,2,4,4,6,6,8,8}; // this routine uses Orlandini-Li et, 2022's  convention
  //std::vector<double> ddp2 = {0,1,3,3,5,5,7,7,1};
  //int conv_type=2;
  // ESRI terra covension:
  //   //     # x-1 x 	x+1
  //   // y-1 ## 32	64	128
  //   // y   ## 16	0	1
  //   // y+1 ## 8	4	2
  //   
  int conv_type=1;
  std::vector<double> ddp1 = {0,64,64,1,1,4,4,16,16};// this routine ESRI  convention
  std::vector<double> ddp2 = {0,32,128,128,2,2,8,8,32};

  //  Li at al, 2022's convection 
  //   //     # x-1 x 	x+1
  //   // y-1 ## 1	2	3
  //   // y   ## 8	0	4
  //   // y+1 ## 7	6	5

  std::vector<double> sigma = {0,1,-1,1,-1,1,-1,1,-1};
  int nncell=9;

  slope_direction(&e[0],nx,ny,&sr[0],&sm[0],&sfacet[0],
                  &tdc[0],&tdd[0],L,ddp1,ddp2,sigma,nncell,conv_type,use_lad);
  bool ok = transverse_deviation(&e[0],&tdc[0],&tdd[0],&sr[0],&sm[0],&sfacet[0],nx,ny,L,

                       &atdc[0], &atdd[0],&atdplus[0],&atdplus0[0], &pflow[0],&has_upstream[0],&kupdate[0],&npids[0],lambda,ddp1,ddp2,sigma,nncell,conv_type,use_lad,max_iters);  
  return ok;
}


SpatRaster  SpatRaster::d8ltd(double lambda,int use_lad,int max_iters,SpatOptions &opt) {
    // DA TESTARE
    SpatRaster out=geometry();
    //std::vector<std::string> oname="watershed";
    //out.setNames(oname);
    int nx=ncol();
    int ny=nrow();
    double Lx=xres();
    //double Ly=yres();
    double L=Lx; // to check 

    //Rprintf("nx=%d ny=%d\n",nx,ny);
    //Rcpp::IntegerVector pOut(nx*ny);
    // https://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027/C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm 
    std::vector<double> e=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
    std::vector<double> pflow(nx*ny,0);
    if (!d8ltd_computation(&e[0],nx,ny,L,lambda,use_lad,max_iters,&pflow[0])) {
      out.addWarning("exceeded the maximum number of iterations in d8ltd/d8lad flow directions computation");
    }
    if (!out.writeStart(opt,filenames())) {
      readStop();
      return out;
    }    
    // convert pflow ....
    // out.writeValues(pOutv,0,ny,0,nx); UNTIL 20220725
    out.writeValues(pflow,0,ny); //,0,nx); // LOOK AT writeValuesGDAL
    out.writeStop();
    // DEALLOC 
    return out;
}  


///// END FLOWDIRECTION


/////// PITFILLER 

// see Grimaldi, S., Nardi, F., Di Benedetto, F., Istanbulluoglu, E., & Bras, R. L. (2007). 
// "A physically-based method for removing pits in digital elevation models." Advances in Water Resources, 30(10), 2151-2158. doi:10.1016/j.advwatres.2006.11.016

/* nx,ny number of cols,rows 
 * pitf pit raster maps
 * npits number of pits
 * e elevaion maps
 * eout modified elevation maps

 * e1 e2 e3  y-1
 * e8 e0 e4  y
 * e7 e6 e5  y+1 */

double calculeZ_pem4ppt(double U, double delta_l, double dx, double D, double z_avg, double beta, double A, double theta, double z_d) {
    double dx2 = dx * dx; // Calcolo di dx al quadrato 
    double numerator = U * delta_l * dx2 + 4 * D * delta_l * z_avg + beta * pow(A, theta) * dx2 * z_d;
    double denominator = beta * pow(A, theta) * dx2 + 4 * D * delta_l;
   // printf("numerator = %f denominator = %f \n",numerator,denominator); 
    double z = numerator / denominator;
    return z;
} 


void pitfiller(int nx,int ny,double ipit, double *pitf,double *e,double *eout,double *flowdir,double L,double *flowacc,
               double U,double D,double beta,double theta_exp) {

  std::vector<int> jv(9,0); /// = {0,1,2,3,4,5,6,7,8};
  std::vector<int> nextv(9,0); /// = {0,1,2,3,4,5,6,7,8};
  std::vector<double> ev(9, 0.0);
  std::vector<double> rlen(9,0); // distances among pixels rescaled with L

  for (int i=0;i<nx*ny;i++) {
    *(eout+i)=*(e+i);
  }
  // loop to detect which pixel of the pit is abjacent with the lowest elevation neighboring pixel
  for (int i=0;i<nx*ny;i++) if (*(pitf+i)==ipit) {
    int y=getRow(nx,ny,i); 
    int x=getCol(nx,ny,i);
    jv[0]=offset(nx,ny,x,y);
    jv[1]=offset(nx,ny,x-1,y-1); //diagonal
    jv[2]=offset(nx,ny,x,y-1);
    jv[3]=offset(nx,ny,x+1,y-1); //diagonal
    jv[4]=offset(nx,ny,x+1,y);
    jv[5]=offset(nx,ny,x+1,y+1); //diagonal
    jv[6]=offset(nx,ny,x,y+1);
    jv[7]=offset(nx,ny,x-1,y+1); //diagonal
    jv[8]=offset(nx,ny,x-1,y); 

    rlen[0]=0;
    rlen[1]=sqrt(2); //diagonal
    rlen[2]=1;
    rlen[3]=sqrt(2); //diagonal
    rlen[4]=1;
    rlen[5]=sqrt(2); //diagonal
    rlen[6]=1;
    rlen[7]=sqrt(2); //diagonal
    rlen[8]=1; 

    //ev[0]=*(e+jv[0]);
    int is_one_pixel=1;
    ev[0] = *(e+jv[0]);
    //printf("ev0=%f _",ev[0]);
    double evavg=0;
    double evdown=*(e+jv[1]);;
    for (int is = 1;is<(int)jv.size();is++) {
      //int yr=getRow(nx,ny,jv[is]); 
      //int xr=getCol(nx,ny,jv[is]);
      //nextv[is]=nextcell_point_conv1(nx,ny,xr,yr,*(flowdir+jv[is]),conv_type);
      // CORRECTION       
      ev[is] = *(e+jv[is]);
 ////     if (is_one_pixel==1) Rprintf("ev[%d]=%f ev[0]=%f evdown=%f\n",is,ev[is],ev[0],evdown);
      evavg=evavg+ev[is];
      if ((ev[is]<=evdown) || (is==1)) evdown=ev[is];
      if (ev[is]<=ev[0]) is_one_pixel=0; // in case there is a pixel lower than the pit pixel!
      if (std::isnan(ev[is])) is_one_pixel=0; // in ase pit is on the boundary and close to nan values. 
    }  

    evavg=evavg/(ev.size()-1);
    // is one pixel?    
    if (is_one_pixel==1) {
    //  double evavg=std::accumulate(ev.begin() + 1, ev.end(),0.0)/(ev.size()-1);
   ////   double evdown=std::min_element(ev.begin()+1, ev.end());
      int isdown=0;
      for (int is = 1;is<(int)jv.size();is++) {
        if (ev[is]==evdown) isdown=is;
      }
      double dx=L;
      double dl=L*rlen[isdown];
      double A_flow=*(flowacc+jv[0])*L*L; // Total Contributing Area [m]
   ////   ev[0]=calculeZ_pem4ppt(double U, double dl, double dx, double D, double z_avg, double beta, double A_flow, double theta_exp, double evdown);
      ev[0]=calculeZ_pem4ppt(U,dl,dx,D, evavg,beta,A_flow,theta_exp,evdown);
    }
    //printf("evout0=%f _\n",ev[0]);
    *(eout+jv[0])=ev[0];
   //20241218 if (is_one_pixel==1) printf("pit number %f : e=%f eout=%f one_pixel=%d\n",ipit,*(e+jv[0]),*(eout+jv[0]),is_one_pixel);
  }  
}


// returns false if the maximum number of iterations was exceeded in the flow directions computation
bool pitfiller_all(int nx,int ny, double *pitf,double *pitftemp,double *e,double *eout,double *flowdir,int niter, 
		double lambda,int use_lad,int max_iters,double L, double U,double D,double beta,double theta_exp) {

  bool ok = true;
  double ipit=-1;
  std::vector<double> nidp_value(nx*ny,0);
  std::vector<int> pnext(nx*ny,0);
  std::vector<double> flowaccu_value(nx*ny,0);
  int pits_on_boundary=0; // pits are considered to be internal the terrain area (dtm)

  for (int i=0;i<nx*ny;i++) {    
    *(eout+i)=*(e+i);
    *(pitftemp+i)=*(pitf+i);
  }

  for (int ir=0;ir<niter;ir++) {
    // Flow Accumlation 
    NextCell(&flowdir[0],nx,ny,&pnext[0]);
    NIDP(&pnext[0],nx,ny,&nidp_value[0]); 
    FlowAccu(&pnext[0],nx,ny,&nidp_value[0],&flowaccu_value[0]);    
    // END flow Accumulation
    for (int i=0;i<nx*ny;i++) {
      for (int i=0;i<nx*ny;i++) {
       // if (*(e+i)!=*(eout+i)) printf("update e[i=%d]=%f eout[i=%d]=%f \n",i,*(e+i),i,*(eout+i));
        *(e+i)=*(eout+i);
      }      
      if (*(pitftemp+i)>0) {

        ipit=*(pitftemp+i);
        pitfiller(nx,ny,ipit,&pitftemp[0],&e[0],&eout[0],&flowdir[0],L,&flowaccu_value[0], U,D,beta,theta_exp); 
      } 
    //20241218  printf("END iteration %07d on %07d \n",ir,niter);
    }
    if (!d8ltd_computation(&eout[0],nx,ny,L,lambda,use_lad,max_iters,&flowdir[0])) { // flow direction computation 
      ok = false;
    }
    pitfinder(&flowdir[0],nx,ny,&pitftemp[0],pits_on_boundary);    ; // pit finder 
   }
  return ok;
}


SpatRaster  SpatRaster::pitfillerm(SpatRaster pits,SpatRaster flowdirs,int niter, double lambda,int use_lad,int max_iters,
  double U,double D,double beta,double theta_exp, // see // see reference doi:10.1016/j.advwatres.2006.11.016)    
  SpatOptions &opt) {
    // DA TESTARE
    SpatRaster out=geometry();
    //std::vector<std::string> oname="watershed";
    //out.setNames(oname);
    int nx=ncol();
    int ny=nrow();
    double Lx=xres();
    double Ly=yres();
    double L=(Lx+Ly)/2; // Lx and Ly are recommended to be equal for a correct usage of the function. 
    //double Ly=yres();
    //double L=Lx; // to check 
  ////  printf("pitfiller\n");
    std::vector<double> e=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
    std::vector<double> pitf=pits.getValues(0,opt);
    std::vector<double> flowdirf=flowdirs.getValues(0,opt);

    std::vector<double> eout(nx*ny,0);
    std::vector<double> pitftemp(nx*ny,0);

    if (!pitfiller_all(nx,ny,&pitf[0],&pitftemp[0],&e[0],&eout[0],&flowdirf[0],niter,lambda,use_lad,max_iters,L,
                  U,D,beta,theta_exp)) { // see reference doi:10.1016/j.advwatres.2006.11.016
      out.addWarning("exceeded the maximum number of iterations in d8ltd/d8lad flow directions computation");
    } 

    if (!out.writeStart(opt,filenames())) {
      readStop();
      return out;
    }

    out.writeValues(eout,0,ny); 
    out.writeStop();
    return out;
}


//// END PITFILLER 

