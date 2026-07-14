#include "spatRaster.h"



#include "catchments_pitfiller.h"
#include "catchments_flow_d8ltd.h"
#include "catchments.h"


/*
 * 
 * nx,ny number of cols,rows 
 * pitf pit raster maps
 * npits number of pits
 * e elevaion maps
 * eout modified elevation maps
 * 
 * 
 */

/*
 * e1 e2 e3  y-1
 * e8 e0 e4  y
 * e7 e6 e5  y+1
 * 
 */

void pitfiller(int nx,int ny,double ipit, double *pitf,double *e,double *eout,double *flowdir,double L,double *flowacc,
               double U,double D,double beta,double theta_exp // see reference doi:10.1016/j.advwatres.2006.11.016
               ) {
  
  std::vector<int> jv(9,0); /// = {0,1,2,3,4,5,6,7,8};
  std::vector<int> nextv(9,0); /// = {0,1,2,3,4,5,6,7,8};
  std::vector<double> ev(9, 0.0);
  std::vector<double> rlen(9,0); // distances among pixels rescaled with L
  double ev1=NAN_PITFILLER;
  double ev2=NAN_PITFILLER;
  double evm=ev1;
  int conv_type=1; // this variable is always 1 see nextcellpoint_conv1
  int jk=WINIT_JK;
  int iss=0;
  
  
  for (int i=0;i<nx*ny;i++) {
    
    *(eout+i)=*(e+i);
    
  }
  // loop to detect which pixel of the pit is abjacent with the lowest elevation neighboring pixel
  /*while (jk==WINIT_JK) */
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
    for (int is = 1;is<jv.size();is++) {
      //int yr=getRow(nx,ny,jv[is]); 
      //int xr=getCol(nx,ny,jv[is]);
      //nextv[is]=nextcell_point_conv1(nx,ny,xr,yr,*(flowdir+jv[is]),conv_type);
      // CORRECTION 
      
      ev[is] = *(e+jv[is]);
 ////     if (is_one_pixel==1) printf("ev[%d]=%f ev[0]=%f evdown=%f\n",is,ev[is],ev[0],evdown);
      
      evavg=evavg+ev[is];
      if (ev[is]<=evdown | is==1) evdown=ev[is];
      
     
      if (ev[is]<=ev[0]) is_one_pixel=0; // in case there is a pixel lower than the pit pixel!
      if (std::isnan(ev[is])) is_one_pixel=0; // in ase pit is on the boundary and close to nan values. 
      
      
    }  
    
    evavg=evavg/(ev.size()-1);
    // is one pixel?
    
    if (is_one_pixel==1) {
    //  double evavg=std::accumulate(ev.begin() + 1, ev.end(),0.0)/(ev.size()-1);
   ////   double evdown=std::min_element(ev.begin()+1, ev.end());
      int isdown=0;
      for (int is = 1;is<jv.size();is++) {
      
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
  
        
void pitfiller_all(int nx,int ny, double *pitf,double *pitftemp,double *e,double *eout,double *flowdir,int niter,double lambda,int max_iters,int use_lad,double L,
                   double U,double D,double beta,double theta_exp) // see // see reference doi:10.1016/j.advwatres.2006.11.016)    
  
                   {
  
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
      
      //20241218 printf("START iteration %07d on %07d \n",ir,niter);
      if (*(pitftemp+i)>0) {
   
        ipit=*(pitftemp+i);
    //  printf("ipit=%f",ipit);
        pitfiller(nx,ny,ipit,&pitftemp[0],&e[0],&eout[0],&flowdir[0],L,&flowaccu_value[0],
        U,D,beta,theta_exp); // see reference doi:10.1016/j.advwatres.2006.11.016
        
        
        
        
        //        for (int b=0;b<nx*ny;b++) {
//         // *(e+b)=*(eout+b);
//        if (*(pitf+b)==ipit) *(pitftemp+b)=0;
      } 
    //20241218  printf("END iteration %07d on %07d \n",ir,niter);
    }
    
    
    
    
    
    d8ltd_computation(&eout[0],nx,ny,L,lambda,use_lad,max_iters,&flowdir[0]); // flow direction computation 
    pitfinder(&flowdir[0],nx,ny,&pitftemp[0],pits_on_boundary);    ; // pit finder 
     
      
      
      
   }

      
      
      
    
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
    double L=Lx; // to check 
    
  ////  printf("pitfiller\n");
    std::vector<double> e=getValues(0,opt); //EC 20211203 //see https://www.delftstack.com/howto/cpp/how-to-convert-vector-to-array-in-cpp/
    std::vector<double> pitf=pits.getValues(0,opt);
    std::vector<double> flowdirf=flowdirs.getValues(0,opt);
    
    std::vector<double> eout(nx*ny,0);
    std::vector<double> pitftemp(nx*ny,0);
    
    pitfiller_all(nx,ny,&pitf[0],&pitftemp[0],&e[0],&eout[0],&flowdirf[0],niter,lambda,use_lad,max_iters,L,
                  U,D,beta,theta_exp); // see reference doi:10.1016/j.advwatres.2006.11.016
    
    
    if (!out.writeStart(opt,filenames())) {
      readStop();
      return out;
    }
    
    out.writeValues(eout,0,ny); 
    out.writeStop();
    return out;
    
    
  }

  

  
  // 
  // 
  // 
  // The PEM4PIT model, which stands for Physical Erosion Model for PIT and flat areas correction, is designed to address issues in Digital Elevation Models (DEMs) related to pits and flat areas. These issues can significantly affect hydrologic modeling and analysis. The model uses a simplified physical representation of erosion processes at the basin scale to correct these inaccuracies.
  // 
  // - **Automatic parameter estimation**: It automatically calibrates parameters to improve the accuracy of DEMs.
  // - **Hydrologic applications**: It enhances the representation of flow directions and contributing areas, which are crucial for hydrologic simulations.
  // - **Erosion processes**: The model accounts for both fluvial and diffusive erosion processes, providing a more realistic depiction of terrain.
  // 
  // Model is based on the steady sediment balance equation:
  //   
  //   \[ 
  // 0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) + \frac{4D}{\Delta x^2} \left( \overset{\smile}{z} - z \right) \tag{1}
  // \]
  // 
  // where:
  //   
  //   - **\(U\)**: uplift
  // - **\(\beta\)**: A coefficient that represents the influence of the erosion process on the terrain.
  // - **\(A\)**: The area contributing to the flow, which affects the erosion rate.
  // - **\(\theta\)**: An exponent that modifies the effect of the contributing area \(A\) on the erosion process.
  // - **\(z\)**: The elevation at a specific point in the DEM.
  // - **\(z_d\)**: The elevation of a downstream point, used to calculate the slope.
  // - **\(\Delta l\)**: The distance between the points \(z\) and \(z_d\), representing the length over which the slope is calculated.
  // - **\(D\)**: A diffusion coefficient that represents the rate of terrain smoothing due to diffusive erosion processes.
  // - **\(\Delta x\)**: The grid spacing in the DEM, used to normalize the diffusion term.
  // - **\(\overset{\smile}{z}\)**: The average elevation of neighboring points, used to calculate the diffusive erosion effect.
  // 
  // ### Solving the Equation for \(z\)
  // 
  // Starting with the equation:
  //   
  //   \[ 
  // 0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) + \frac{4D}{\Delta x^2} \left( \overset{\smile}{z} - z \right) 
  //   \]
  // 
  // 1. Rearrange the terms to isolate \(z\):
  //   
  //   \[ 
  // \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) = U + \frac{4D}{\Delta x^2} \left( \overset{\smile}{z} - z \right) 
  //   \]
  // 
  // 2. Distribute \(\beta A^{\theta}\) and \(\frac{4D}{\Delta x^2}\):
  //   
  //   \[ 
  // \beta A^{\theta} \frac{z}{\Delta l} - \beta A^{\theta} \frac{z_d}{\Delta l} = U + \frac{4D \overset{\smile}{z}}{\Delta x^2} - \frac{4D z}{\Delta x^2} 
  // \]
  // 
  // 3. Combine like terms involving \(z\):
  //   
  //   \[ 
  // \beta A^{\theta} \frac{z}{\Delta l} + \frac{4D z}{\Delta x^2} = U + \frac{4D \overset{\smile}{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l} 
  // \]
  // 
  // 4. Factor out \(z\):
  //   
  //   \[ 
  // z \left( \beta A^{\theta} \frac{1}{\Delta l} + \frac{4D}{\Delta x^2} \right) = U + \frac{4D \overset{\smile}{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l} 
  // \]
  // 
  // 5. Solve for \(z\):
  //   
  //   \[ 
  // z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \overset{\smile}{z} + \beta A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta l} 
  // \]
  // 
  // So, the solution for \(z\) is:
  //   
  //   \[ 
  // z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \overset{\smile}{z} + \beta A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta l} 
  // \]
  // 
  // 
  // 
  // see reference: 
  // Grimaldi, S., Nardi, F., Di Benedetto, F., Istanbulluoglu, E., & Bras, R. L. (2007). 
  // "A physically-based method for removing pits in digital elevation models."
  // Advances in Water Resources, 30(10), 2151-2158. doi:10.1016/j.advwatres.2006.11.016
  //
  //
  
  double calculeZ_pem4ppt(double U, double delta_l, double dx, double D, double z_avg, double beta, double A, double theta, double z_d) {
    double dx2 = dx * dx; // Calcolo di dx al quadrato
    
    // Print input variables
    // printf("U = %f\n", U);
    // printf("delta_l = %f\n", delta_l);
    // printf("dx = %f\n", dx);
    // printf("D = %f\n", D);
    // printf("z_avg = %f\n", z_avg);
    // printf("beta = %f\n", beta);
    // printf("A = %f\n", A);
    // printf("theta = %f\n", theta);
    // printf("z_d = %f\n", z_d);
    // printf("dx^2 = %f\n", dx2);
    
    
    double numerator = U * delta_l * dx2 + 4 * D * delta_l * z_avg + beta * pow(A, theta) * dx2 * z_d;
    double denominator = beta * pow(A, theta) * dx2 + 4 * D * delta_l;
   // printf("numerator = %f denominator = %f \n",numerator,denominator);
    
    
    double z = numerator / denominator;
    return z;
  } 
    
  
  
 











