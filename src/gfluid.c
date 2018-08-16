#include "config.h"
#include "gfluid.h"

control_params 	Control;
evolve_params 	EvolveConfig;
sim_data 	DataCurr, DataPrev;
sim_data 	DataSpectrum, DataSurface;


/* main function */
int main (int argc, char **argv) {
  Control.DataPtrCurr = &DataCurr;
  Control.DataPtrPrev = &DataPrev;
  Control.EvolvePtr = &EvolveConfig;
  
  evaluate_jacobi_sn();

  return 0;
}
