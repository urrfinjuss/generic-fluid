#define MPFR_WANT_FLOAT128

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <quadmath.h>
#include <complex.h>
#include <fftw3.h>
#include <mpfr.h>
#include <mpc.h>

//#define PIq 	acosq(-1.0Q)
#define PIq	M_PIq
#define PIl	acosl(-1.0L)
#define MPFR_PREC	512

typedef __float128		long_double_t;
typedef fftwq_complex		long_complex_t;
typedef struct data_array	*data_ptr;
typedef struct conformal_map	*cmap_ptr;

#include "headers/mapping.h"
#include "headers/io.h"

typedef struct stepping_parameters {
  unsigned long nsteps, cur_step, dmp_cnt;
  __float128 	cfl, max_cfl, dt;
  __float128	final_time;
} evolve_params, *evolve_params_ptr;


/*  Note: Move to memory.h or array.h  */

typedef struct data_array {
  unsigned long		N;
  long_complex_t	*Q, *V;
  long_complex_t	*Z, *Phi;
  //long_complex_t	z0, beta; // what is it and why?
  //long_double_t		*q;
  cmap_ptr		map;
  //long_complex_t	*u, *du;
  //long_double_t		q0, u0, l;
  long_double_t		Volume, Hamiltonian, SurfaceEnergy;
  __float128 		time;
} sim_data, *data_ptr;

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

#include "headers/memory.h"

/*  END  */



/* Note to self: remember top to bottom this time */
#include "headers/special_functions.h"
#include "headers/math.h"
#include "headers/timemarching.h"

#include "headers/debug.h"
//#include "headers/array_func.h"


#include "headers/messages.h"

/* Control Parameters: move to parameters.h  */

typedef struct control_parameters_array{
  data_ptr		DataPtrCurr, DataPtrPrev;
  evolve_params_ptr 	EvolvePtr;
  __float128		Sigma;
  __float128		Gravity;
  __float128		transform_kc;
  char 			run_name[80];
  char			res_name[80];
  char			data_path[80];
} control_params, *control_params_ptr;

/* Global Variable: here is fine */ 
extern sim_data		DataCurr, DataPrev;
extern sim_data         DataSpectrum, DataSurface;
extern control_params	Control;
extern evolve_params	EvolveConfig;
extern cmap this_map;
extern cmap last_map;

