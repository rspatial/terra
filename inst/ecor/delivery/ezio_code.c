/*
  * Function:   resizeQueue
  * Scope:      resize (doubling) an existing queue to store integer offsets of raster cells
  * Parameters: int* q pointer to existing queue
  *             int n  current number of elements
  * Return a pointer to the resized queue, while preserving previous values
  */
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
  
  /*
   * Scope: compute watershed upstream of point i,j
   * p: pointer to an integer array storing a 2D raster
   * nx,ny: number of cells in x (longitude) and in y (latitude)
   * x,y: indexes of cell upstream of which the watershed must be computed
   * NOTE: the function is an iterative version of the previous recursive version
   * Current version implements a dinamic queue, up to the basin boundaries
   */
  void watershed_v2_rq(int* p, int nx, int ny, int x, int y, int* pOut)
  {
    int* q;           // A pointer to a queue of raster cells (offset in memory) to be processed
    int qSize = 50;   // Starting queue size, that can be dinamically incremented if needed
    int delta;        // Offset in memory from base queue address of a raster cell
    int n = 0;        // Number of raster cells to be processed in queue
    int nLoop = 0;    // Counter for loops over cells
    
    printf("DEBUG: col=%d,row=%d\n", x, y);
    
    q = (int*)CPLMalloc(sizeof(int)*qSize);
    
    // printf("TEST Row, cell 60: %d,%d\n", getRow(nx,ny,60), getCol(nx,ny,60));
    // Set raster cell in the output file
    delta = offset(nx, ny, x, y);
    *(pOut + delta) = 1;
    // Store raster cell offset in the queue  and update number of elements
    q[0] = delta;
    n++;
    
    printf("BEFORE n=%d and size(n)=%d\n", n, sizeof(n));
    
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
        puts("Press any key to continue ... ");
        getchar();
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
  