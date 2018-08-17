#include "config.h"
#include "gfluid.h"

control_params 	Control;
evolve_params 	EvolveConfig;
sim_data 	DataCurr, DataPrev;
sim_data 	DataSpectrum, DataSurface;
cmap		this_map;
cmap		last_map;

/* main function */
int main (int argc, char **argv) {
  Control.DataPtrCurr = &DataCurr;
  Control.DataPtrPrev = &DataPrev;
  Control.EvolvePtr = &EvolveConfig;
  
  ffluid_read_cl_arguments(argc, argv);
  
  set_HTmap(&this_map);
  gfluid_write_map("ht.cmg");
  
  //elliptic_demo();

  return 0;
}
