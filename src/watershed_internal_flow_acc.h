// C/C++ code
// Author: Emanuele Cordano
// Date: February 2021
// Flow accumulation C routines
//

void NextCell(double* p, int nx, int ny,int* pnext);

void NIDP(int* pnext, int nx, int ny,int* nidp_value); 
// FlowAccu algorithm 5 
// Reference: https://link.springer.com/article/10.1007/s11707-018-0725-9
void FlowAccu(int* pnext, int nx, int ny,double* nidp_value,double* flowaccu_value);

void FlowAccu_weight(int* pnext, int nx, int ny,double* nidp_value,double* flowaccu_value,double* weight);



