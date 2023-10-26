
//
// Scope: compute watershed upstream of point i,j
// p: pointer to an integer array storing a 2D raster
// nx,ny: number of cells in x (longitude) and in y (latitude)
// x,y: indexes of cell upstream of which the watershed must be computed
// NOTE: the function is an iterative version of the previous recursive version
// Current version implements a dinamic queue, up to the basin boundaries
//
// void watershed_v2(int* p, int nx, int ny, int x, int y, int* pOut)
void pitfinder(double* p, int nx, int ny, double* pOut)
{
  int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
//  int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
  int delta;        // Offset in memory from base queue address of a raster cell
  int n = 0;        // Number of raster cells to be processed in queue
//  int nLoop = 0;    // Counter for loops over cells
  int x,y;
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
  q[0] = delta;
  n++;
  
  printf("BEFORE n=%d and size(n)=%d\n", n, sizeof(n));
  for (int i = 0; i < nx*ny; i++) {
    
    *(pOut+i)=0;  
    // ## 32	64	128
    // ## 16	x	1
    // ## 8	4	2
    
    
    // ## 2 4 8
    // ## 1 y 16
    // ##128   64   32    
    
    x = getCol(nx, ny,i);   // ATTENTION: base 0 or 1?
    y = getRow(nx, ny,i);
    
    if (*(p+i) == 1) {
      
      if (inRaster(nx, ny, x + 1, y) && *(p+i+offset(nx, ny, x + 1, y)) == 16) {
        
        *(pOut+i)=1;
      } 
    } else if (*(p+i) == 2) {
    
      if (inRaster(nx, ny, x + 1, y-1) && *(p+i+offset(nx, ny, x + 1, y-1) == 32) {
      
      *(pOut+i)=1;
      } 
      
    } else if (*(p+i) == 4) {
      
      if (inRaster(nx, ny, x, y-1) && *(p+i+offset(nx, ny, x, y-1) == 64) {
        
        *(pOut+i)=1;
      } 
            
    } else if (*(p+i) == 8) {
      
      if (inRaster(nx, ny, x, y-1) && *(p+i+offset(nx, ny, x - 1, y-1) == 128) {
        
        *(pOut+i)=1;
      } 
    } else if (*(p+i) == 16) {
      
      if (inRaster(nx, ny, x - 1, y-1) && *(p+i+offset(nx, ny, x - 1, y) == 1) {
        
        *(pOut+i)=1;
      } 
    } else if (*(p+i) == 32) {
      if (inRaster(nx, ny, x + 1, y-1) && *(p+i+offset(nx, ny, x - 1, y+1) == 2) {
        
        *(pOut+i)=1;
      } 
    } else if (*(p+i) == 64) {
      if (inRaster(nx, ny, x + 1, y-1) && *(p+i+offset(nx, ny, x - 1, y+1) == 2) {
        
        *(pOut+i)=1;
      } 
    } else if (*(p+i) == 128) {
      if (inRaster(nx, ny, x + 1, y-1) && *(p+i+offset(nx, ny, x + 1, y+1) == 2) {
        
        *(pOut+i)=1;
      } 
    }
    
    // SI CONTINUA DA QUI 
    
    
    
  //   
  //   if (inRaster(nx, ny, x + 1, y) && *(p+i+offset(nx, ny, x + 1, y)) == 16) {
  //     
  //     if (*(p+i) == 1) { 
  //       *(pOut+i)=1;
  //       }
  //   }
  //   
  //   // and so omn
  //   if (*(p+i)== pdown=*(&p[i] + offset(nx, ny, x + 1, y)))
  //   
  //   
  //   if (inRaster(nx, ny, x + 1, y) && *(p + offset(nx, ny, x + 1, y)) == 16) {
  //     
  //     
  //   }
  //   
  //   
  //   
  // }
  // 
  // q[i] = q[i + 1];
  // // Process all pending (until any) raster cells in the queue
  // while (n > 0) {
  //   //printf("DEBUG: IN THE LOOP n=%d\n", n);
  //   nLoop++;
  //   if (nLoop % 100000 == 0) printf("%d\n", nLoop);  // Print number of internal loops
  //   
  //   
  //   
  //   //printf("%d\n", n);
  //   // Pick up top raster cell
  //   x = getCol(nx, ny, q[0]);   // ATTENTION: base 0 or 1?
  //   y = getRow(nx, ny, q[0]);
  //   //printf("col=%d row=%d\n", x, y);
  //   if (inRaster(nx, ny, x + 1, y) && *(p + offset(nx, ny, x + 1, y)) == 16) {
  // 
  //   
  //   // Queue is just full, only 10 raster cells to accomodate
  //   if (n > (qSize-10)) {
  //     q = resizeQueue(q, qSize);
  //     qSize *= 2;
  //     // puts("Press any key to continue ... ");
  //     // getchar();
  //   }
  //   
  //   // Investigate D8 raster cells all around the cell under investigation
  //   // Bounding raster cell located to the E
  //   if (inRaster(nx, ny, x + 1, y) && *(p + offset(nx, ny, x + 1, y)) == 16) {
  //     delta = offset(nx, ny, x + 1, y);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the SE
  //   if (inRaster(nx, ny, x + 1, y + 1) && *(p + offset(nx, ny, x + 1, y + 1)) == 32) {
  //     delta = offset(nx, ny, x + 1, y + 1);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the S
  //   if (inRaster(nx, ny, x, y + 1) && *(p + offset(nx, ny, x, y + 1)) == 64) {
  //     delta = offset(nx, ny, x, y + 1);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the SW
  //   if (inRaster(nx, ny, x - 1, y + 1) && *(p + offset(nx, ny, x - 1, y + 1)) == 128) {
  //     delta = offset(nx, ny, x - 1, y + 1);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the W
  //   if (inRaster(nx, ny, x - 1, y) && *(p + offset(nx, ny, x - 1, y)) == 1) {
  //     delta = offset(nx, ny, x - 1, y);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the NW
  //   if (inRaster(nx, ny, x - 1, y - 1) && *(p + offset(nx, ny, x - 1, y - 1)) == 2) {
  //     delta = offset(nx, ny, x - 1, y - 1);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the N
  //   if (inRaster(nx, ny, x, y - 1) && *(p + offset(nx, ny, x, y - 1)) == 4) {
  //     delta = offset(nx, ny, x, y - 1);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   // Bounding raster cell located to the NE
  //   if (inRaster(nx, ny, x + 1, y - 1) && *(p + offset(nx, ny, x + 1, y - 1)) == 8) {
  //     delta = offset(nx, ny, x + 1, y - 1);
  //     *(pOut + delta) = 1;
  //     q[n] = delta;
  //     n++;
  //   }
  //   
  //   //printf("DEBUG AT THE END - n=%d\n",n);
  //   
  //   // Just completed the analysis on current raster cell.
  //   // Raster cells in queue are shifted up of one position, current raster cell being
  //   // the facto removed from the queue, and the number of elements is updated
  //   for (int i = 0; i < n; i++) q[i] = q[i + 1];
  //   n--;
  // }
  // 
  // CPLFree(q);  // Frees memory allocated to the queue (whether is the original or a resized one)
}



// 
//   
// //
// //
// int flowdir(int i, double* p, int nx, int ny)
// {
//   double val=*(p+i);
//   // ## 32	64	128
//   // ## 16	val	1
//   // ## 8	4	2
//   int out;
//   
//   
//   if (val==1) 
//   
//   
//   y * nx + x; //according to original Ezio's code BY ROWS 
//   //return x * ny + y; // according offset defintion in Rccp for IntergerMatrix // BY COLS
// }
//   
// 
// 
// 
// 
// 
// 
// 
