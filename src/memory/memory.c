#include "gfluid.h"


void ffluid_memory_module() {
  printf("Module ffluid/memory.h:\n");
  printf("memory.c:\tffluid_memory_module ffluid_init_data\n");
}

void ffluid_init_data(data_ptr in) {
  in->R = fftwq_malloc(in->N*sizeof(long_complex_t));
  in->V = fftwq_malloc(in->N*sizeof(long_complex_t));
  in->u = fftwq_malloc(in->N*sizeof(long_complex_t));
  in->du = fftwq_malloc(in->N*sizeof(long_complex_t));
  in->q = fftwq_malloc(in->N*sizeof(long_double_t));
}

void ffluid_alloc_aux_array(aux_data_ptr in, unsigned long NArrays, unsigned long NElements) {
  in->NArrays = NArrays;
  in->NElements = NElements;
  printf("\tAllocate %lu arrays of %lu elements:", NArrays, NElements);
  in->X = fftwq_malloc(NArrays*sizeof(long_complex_t *));
  in->Y = fftwq_malloc(NArrays*sizeof(long_complex_t *));
  for (unsigned long j = 0; j < NArrays; j++) {
    in->X[j] = fftwq_malloc(NElements*sizeof(long_complex_t));
    in->Y[j] = fftwq_malloc(NElements*sizeof(long_complex_t));
  }
  printf(" ... memory allocated\n");
}

void ffluid_dealloc_aux_array(aux_data_ptr in) {
  for (unsigned long j = 0; j < in->NArrays; j++) {
    fftwq_free(in->X[j]);
    fftwq_free(in->Y[j]);
  }
  fftwq_free(in->X);
  fftwq_free(in->Y);
}

void ffluid_alloc_fft_plans(aux_data_ptr in, fft_list_ptr out) {
  out->NFFTs = in->NArrays;
  printf("\tInitializing %lu FFT plans:", in->NArrays);
  out->fp = fftwq_malloc(in->NArrays*sizeof(fftwq_plan));
  out->bp = fftwq_malloc(in->NArrays*sizeof(fftwq_plan));
  for (unsigned long j = 0; j < out->NFFTs; j++) {
    out->fp[j] = fftwq_plan_dft_1d(in->NElements, in->X[j], in->Y[j], FFTW_BACKWARD, FFTW_ESTIMATE);
    out->bp[j] = fftwq_plan_dft_1d(in->NElements, in->Y[j], in->X[j], FFTW_FORWARD, FFTW_ESTIMATE);
  }
  printf(" ... initialized\n");
}


