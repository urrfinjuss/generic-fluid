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
 

  /*  Various inits */
  ffluid_read_cl_arguments(argc, argv);
  ffluid_set_initial_data(&DataCurr);
  ffluid_alloc_equations();

  gfluid_emergency_init(&DataCurr, 1.0E-11IQ, 2.0E-2IQ, 1.0E-1IQ, 0.2Q);
  gfluid_write_map(DataCurr.map, DataCurr.N, "ht.cmg");
  
  gfluid_get_surface(&DataCurr);
  ffluid_write_natural(&DataCurr, "natural_000.txt");
  ffluid_write_raw(&DataCurr, "raw_000.txt");
  
  //elliptic_demo();

  return 0;
}
