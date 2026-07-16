// C/C++ code
// Author: Emanuele Cordano
// Date: February 2021
// Flow direction D8-LTD algorithm implementation (under development)
//


//See this former code: 
//line 217 qui: https://github.com/TheHortonMachine/hortonmachine/blob/master/hmachine/src/main/java/org/hortonmachine/hmachine/modules/geomorphology/draindir/OmsDrainDir.java per la prescisione 🙂

// TO BE IMPLEMENTED

#include "spatRaster.h"
#include "NA.h"
#include <cmath> // per la funzione pow
#include "catchments.h"

#define NODATA_2 -999999

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

/*
 * 
 * nx,ny number of cols,rows 
 * x,y    row,col
 * pv flow direction value 
 * conv_type convention type used , see code
 * 
 
 * 
 */
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
  
  /*
   * 
   * e elevation
   * nx,ny number of cols,rows
   * sr SR in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * sm SM in Li et al,2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * sfacet drainage facet 
   * tdc TD_c in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   * tdd TD_c in Li et, 2022 (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021WR030948)
   *
   * 
   */

void slope_direction(double* e, int nx, int ny, double *sr,double *sm,int *sfacet,
                     double* tdc, double *tdd,double L,
                     std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,int conv_type,int use_lad) {
  
  int i;
  ///int L=1; // dx=dy=L=1
  int ncell=nx*ny;
  double x,y;
 // int efacet=0;
  double e0,e1,e2;
  int facet=0;
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
    double flow_angle_tan=0;
    
    
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
      
      
      
      
      
      
      // nothing 
    //}
    
  ///  slope_mgn_temp=pow(pow(e0-e1,2)+pow(e1-e2,2),0.5)/L;
    //double flow_angle_tan=(e1-e2)/(e0-e1); //ec 20251023
   /// double flow_angle_tan=((e1-e2)/(e0-e1)); // removed abs 20251011  
    //*sigma[facetc] DA VEDERE 20251110
    //*
    //slope_mgn=slope_mgn_temp; //pow(pow(mean_e-e1,2)+pow(e1-e2,2),0.5)/L;
    //if (e0<mean_e_temp) { // uncommented on EC 20250717
     //  
    //   slope_mgn_temp=0;
      
      
    //}  // else // commented on EC 20251110
    // if (flow_angle_tan>1)  {
    //   slope_mgn_temp=((e0-e2)/(L*sqrt(2))); //20251109 if slope direction is outside of the facet  cannot be considerated
    //   ////slope_mgn_temp=0;
    // } else if (flow_angle_tan<0) {
    //   slope_mgn_temp=((e0-e1)/L);  //20251112  not abs!!
    // } // 20251022
    
    //if (slope_mgn_temp<0) slope_mgn_temp=0;
    
    if (slope_mgn_temp>slope_mgn) {
        slope_mgn=slope_mgn_temp;
        mean_e=mean_e_temp;
        facet=j;
        
    } else  if (slope_mgn_temp==slope_mgn) {
       
          //Rprintf("\n i=%d j=%d ",i,j);
          //Rprintf("slope_mgn_temp=%f slope_mgn=%f mean_e_temp=%f mean_e=%f",slope_mgn_temp,slope_mgn,mean_e_temp,mean_e);
          //Rprintf("facet=%d \n",facet); 
         
         if (mean_e_temp<mean_e) {
          slope_mgn=slope_mgn_temp;
          mean_e=mean_e_temp;
          facet=j;
          
          
          } else if (mean_e_temp==mean_e){
              facet=0; // uncommented on 20250109 // more actions to do // 2026 01 18 
          
          }
        
        
        
     }
   //  Rprintf("ab i=%d x=%f y=%f j=%d e0=%f e1=%f e2=%f",i,x,y,j,e0,e1,e2);
   //  Rprintf("slope_mgn_temp=%f slope_mgn=%f mean_e_temp=%f mean_e=%f \n",slope_mgn_temp,slope_mgn,mean_e_temp,mean_e);
  //   Rprintf("facet=%d \n\n",facet); 
     
      
    }
    
    double e_temp=e0;
    if (facet==0)  for (int j=0; j<nncell; j++){
      
      //if (ddp1[j]==1) {
      /// Rprintf("ba2 j=%d facet=%d i=%d e0=%f e1=%f e2=%f\n",j,facet,i,e0,e1,e2);
      e1=*(e+nextcell_point_conv1(nx,ny,x,y,ddp1[j],conv_type));
      e2=*(e+nextcell_point_conv1(nx,ny,x,y,ddp2[j],conv_type));
      mean_e_temp=(e0+e1+e2)/3;
      if (!(std::isnan(e1) | std::isnan(e2))) { // to rewrite 
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
   // Rprintf("j=%d ",j);
   //cc20260617 Rprintf("i=%d x=%f y=%f facet=%d e0=%f e1=%f e2=%f \n",i,x,y,facet,e0,e1,e2); 
   ///ghj double slope_test=(std::atan((e1-e2)/(e0-e1)));
   //cc20260617 Rprintf("slope_mgn=%f   slope_test=%f e0=%f e1=%f e2=%f  \n\n",*(sm+i),*(sr+i),e0,e1,e2);
    
    
    
  
}
      
      
      
      
      
}
 
 
void transverse_deviation(double *e, double *tdc, double *tdd,double *sr,double *sm, int *sfacet,int nx, int ny, double L,
                          
                          double *atdc, double *atdd, double *atdplus,double *atdplus0,
                          double *pflow,int *has_upstream,int *kupdate,double *nidps,
                          double lambda,
                          std::vector<double> ddp1,std::vector<double> ddp2,std::vector<double> sigma,int nncell,int conv_type,int use_lad,int max_iters)    {   
   
 
  int x,y;

  int k=1;
  int niter=nx*ny;
  int exit_cond=0;
  int exit_cond1=0;
  int facet,nextc,nextd;
  int nextp=0;
  double atdplus_temp;
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
   } else if ((abs(*(atdc+i))<abs(*(atdd+i))) | ((abs(*(atdc+i))<=abs(*(atdd+i))) & (use_lad==1))){
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
      //
    
      //
 

    // Rprintf("HERE!!! x=%d y=%d i=%d facet=%d e0=%f e1=%f e2=%f L=%f atdc=%f atdd=%f\n",x,y,i,facet,e0,e1,e2,L,atdplus_nextc,atdplus_nextd);
    // Rprintf("HERE!!! x=%d y=%d i=%d facet=%d facetc=%d facetd=%d e0=%f e1=%f e2=%f L=%f tdc=%f tdd=%f tdcs=%f tdds=%f\n",x,y,i,facet,facetc,facetd,e0,e1,e2,L,*(tdc+nextc),*(tdd+nextd),*(tdc+nextc)*sigma[facetc],*(tdd+nextd)*sigma[facetd]);
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
    } else if ((abs(atdplus_nextpc)<abs(atdplus_nextpd)) | ((abs(atdplus_nextpc)<=abs(atdplus_nextpd)) & (use_lad==1))){
        
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
   //21 Rprintf("zzx i=%d pflow=%f %f\n",i,pflow_estimate,*(pflow+i));
  //  if (pflow_estimate==*(pflow+i)) {
    //  exit_cond=1;
  //21    Rprintf("zzz pflow=%f pfow=%f\n",pflow_estimate,*(pflow+i));
  //21    Rprintf("zzz i=%d nextp=%d su j=%d facet=%d \n",i,nextp,j,facet);
      
 //   }  else {
      
 //     exit_cond=0;
 //   }///else if (nextp!=i) {
    
    // work here 20260122
    // condizione sopre ok se non cambia direzione ??? se e uguale e cambia dierezione???
    if ((abs(atdplus_temp)>=abs(*(atdplus+nextp)))) { // 20260119    if (abs(atdplus_temp)>=abs(*(atdplus+i))){ // 20260119
        // ADD A CONTROL (kupdate+nextp )
      *(atdplus+nextp)=atdplus_temp; // corrected on 20260429
      *(pflow+nextp)=pflow_estimate;
      *(kupdate+nextp)=cnt;
    //  *(kupdate+i)=cnt;
        // aggiorna il pfolw ??
        
    //} 
    // else  if (*(pflow+nextp)!=pflow_estimate) {
    //   int nextp_old=nextcell_point_conv1(nx,ny,x,y,*(pflow+nextp),conv_type);
    //   *(kupdate+nextp_old)=0;
    //   *(pflow+nextp)=pflow_estimate;
    //   *(kupdate+nextp)=cnt;
      /// to do something 
      
      
    } else {
      
      //*(pflow+nextp)=pflow_estimate;
      //*(kupdate+nextp)=cnt; 
      
    };{  //}if (CONDITION) {
     //EC20260617  *(pflow+nextp)=pflow_estimate;
     //EC20260617 *(kupdate+nextp)=cnt;
    //  nextp=nextp0;
    //  exit_cond=3; 
    }
    ///exit_cond=0;
    //*(pflow+i)=pflow_estimate;
   // *(pflow+nextp)=pflow_estimate;
   // *(kupdate+i)=cnt; // corrected on 20260429
  //  *(kupdate+nextp)=cnt;
   //// 20241030 
   //cc 20260617 Rprintf("x=%d y=%d atdd=%f atdc=%f tdd=%f tdc=%f sm=%f sr=%f pflow=%f kup=%d i=%d j=%d nextp=%d\n facet=%d e0=%f e1=%f e2=%f \n",x,y,*(atdd+i),*(atdc+i),*(tdd+i),*(tdc+i),*(sm+i),*(sr+i),*(pflow+i),*(kupdate+i),i,j,nextp,facet,e0,e1,e2);
   //// Rprintf("nextq=%d nextp=%d i=%d\n",nextq,nextp,i);
    if (nextp!=i) { // 20260429
      
      i=nextp;
      nextp=nextq;
      exit_cond=0;
    } else {
      
      exit_cond=2;
    }
      
      
 // }
    
    //}  else {
    //  Rprintf("L2 i=%d nextp=%d su j=%d facet=%f pflow=%f\n",i,nextp,j,facet,*(pflow+i));
  //    exit_cond=1;
  //  }
  ///  Rprintf("j=%d su i=%d exit_cond=%d lambda=%f facet=%d\n",j,i,exit_cond,lambda,facet);
    
  } while (exit_cond==0);
  cnt1++;
  exit_cond1=2;
  
  // for statemant to verify kupdate!=0
  }
  for (int j = 0; j < nx*ny; j++) {
    
    if (*(kupdate+j)==0) exit_cond1=0;
  }
  if (cnt1>max_iters) {
    Rprintf("\nExceeding number of iterations in d8ltd/d8lad flow directions computation\n");
    exit_cond1=2;       
    
  }
  
 } while  (exit_cond1==0);
  // NOVALUE
  for (int i=0;i<nx*ny;i++) {
  //  
  
    double e0=*(e+i);
  
    if (std::isnan(e0)) {
      *(pflow+i)=e0;
   }
  

 
  
  
  
}
 


    
    
    
    

  
  
}


void d8ltd_computation(double *e,int nx,int ny,double L,double lambda,int use_lad,int max_iters,double*pflow) {
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
  
  
  
  
  
  transverse_deviation(&e[0],&tdc[0],&tdd[0],&sr[0],&sm[0],&sfacet[0],nx,ny,L,
                       
                       &atdc[0], &atdd[0],&atdplus[0],&atdplus0[0], &pflow[0],&has_upstream[0],&kupdate[0],&npids[0],lambda,ddp1,ddp2,sigma,nncell,conv_type,use_lad,max_iters);
  
  
  // NOVALUE
  //for (int i=0;i<nx*ny;i++) {
  //  
  
  //  double e0=*(e+i);
  
  //  if (e0==0) {
  //    *(pflow+i)=*(e+i);
  // }
  
}







  
  
SpatRaster  SpatRaster::d8ltd(double lambda,int use_lad,int max_iters,SpatOptions &opt) {
    // DA TESTARE
    SpatRaster out=geometry();
    //std::vector<std::string> oname="watershed";
    //out.setNames(oname);
    int nx=ncol();
    int ny=nrow();
    double Lx=xres();
    double Ly=yres();
    double L=Lx; // to check 

    //Rprintf("nx=%d ny=%d\n",nx,ny);
    //Rcpp::IntegerVector pOut(nx*ny);
    // https://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027/C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm 
    std::vector<double> e=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
    std::vector<double> pflow(nx*ny,0);
    d8ltd_computation(&e[0],nx,ny,L,lambda,use_lad,max_iters,&pflow[0]);
    
  
    
    
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
  

