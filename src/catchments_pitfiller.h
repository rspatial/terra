// C/C++ code
// Author: Emanuele Cordano
// Date: February 2021
// Pit Filler
//

// TO BE IMPLEMENTED

//#define NODATA0 -9999


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

double calculeZ_pem4ppt(double U, double delta_l, double dx, double D, double z_avg, double beta, double A, double theta, double z_d); 
