#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <quadmath.h>
#include <complex.h>
#include <fftw3.h>

#define PIq 	acosq(-1.0Q)
#define PIl	acosl(-1.0L)

typedef __float128	long_double_t;
typedef fftwq_complex	long_complex_t;

typedef struct stepping_parameters {
  unsigned long nsteps, cur_step, dmp_cnt;
  __float128 	cfl, max_cfl, dt;
  __float128	final_time;
} evolve_params, *evolve_params_ptr;

typedef struct data_array {
  unsigned long		N;
  long_complex_t	*R, *V;
  long_complex_t	*u, *du;
  long_complex_t	z0, beta;
  long_double_t		*q;
  long_double_t		q0, u0, l;
  long_double_t		Volume, Hamiltonian, SurfaceEnergy;
  __float128 		time;
} sim_data, *data_ptr;

#include "headers/io.h"

typedef struct aux_array {
  unsigned long		NElements;
  unsigned long		NArrays;
  long_complex_t	**X;
  long_complex_t	**Y;
} aux_data, *aux_data_ptr;

typedef struct fft_array {
  unsigned long		NFFTs;
  fftwq_plan		*fp;
  fftwq_plan		*bp;
} fft_list, *fft_list_ptr;

/* Note to self: remember top to bottom this time */
#include "headers/memory.h"
#include "headers/mapping.h"
#include "headers/special_functions.h"
#include "headers/debug.h"
//#include "headers/array_func.h"
#include "headers/messages.h"

typedef struct control_parameters_array{
  data_ptr		DataPtrCurr, DataPtrPrev;
  evolve_params_ptr 	EvolvePtr;
  long_double_t		Sigma;
  char 			run_name[80];
  char			res_name[80];
  char			data_path[80];
} control_params, *control_params_ptr;

/* Global Variable */ 
extern sim_data		DataCurr, DataPrev;
extern sim_data         DataSpectrum, DataSurface;
extern control_params	Control;
extern evolve_params	EvolveConfig;


