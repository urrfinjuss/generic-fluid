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
  ffluid_alloc_operators();
  ffluid_setup_stepping();

  gfluid_emergency_init_dipole(&DataCurr, 1.0E-1IQ, 1.5E-1IQ, 2.5E+0Q, 0.2Q);
  //gfluid_emergency_init_mode(&DataCurr, 5.0E-1Q, 0.2Q);
  gfluid_write_map(DataCurr.map, DataCurr.N, "ht.cmg");
  
  gfluid_data_init_copy(&DataCurr, &DataSpectrum);
  gfluid_make_spectrum(&DataSpectrum);
  ffluid_write_spectrum(&DataSpectrum, "spectrum_000.txt");
  
  gfluid_get_surface(&DataCurr);
  ffluid_write_natural(&DataCurr, "natural_000.txt");
  ffluid_write_raw(&DataCurr, "raw_000.txt");
  

  

  /*
  cmap_ptr map = DataCurr.map;
  for (long int j = 0; j < DataCurr.N; j++) {
    DataCurr.Q[j] = 1.0Q + 1.0Q*sinq(map->u[j]);
  }
  gfluid_projector(DataCurr.Q, DataCurr.V);
  ffluid_write_raw(&DataCurr, "projector_000.txt");
  */

  ffluid_init_runge_kutta_4();
  ffluid_evolve();
  ffluid_write_raw(&DataCurr, "raw_001.txt");
  
  gfluid_get_surface(&DataCurr);
  ffluid_write_natural(&DataCurr, "natural_001.txt");
  gfluid_data_init_copy(&DataCurr, &DataSpectrum);
  gfluid_make_spectrum(&DataSpectrum);
  ffluid_write_spectrum(&DataSpectrum, "spectrum_001.txt");

  

  return 0;
}
