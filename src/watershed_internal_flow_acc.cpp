//#include <Rcpp.h>
//using namespace Rcpp;
#include "watershed_internal.h"
#include "spatRaster.h"
#include "watershed_internal_flow_acc.h"


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
