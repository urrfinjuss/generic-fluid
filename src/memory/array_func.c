#include "gfluid.h"

void ffluid_arrayfunc_module() {
  printf("Module ffluid/array.h:\n");
  printf("array_func.c:\tffluid_arrayfunc_module ffluid_data_copy ffluid_data_fma\n");
}

void ffluid_data_copy(data_ptr in, data_ptr out) {
  memcpy(out->Q,   in->Q,   in->N*sizeof(long_complex_t));
  memcpy(out->V,   in->V,   in->N*sizeof(long_complex_t));
  memcpy(out->Z,   in->Z,   in->N*sizeof(long_complex_t));
  memcpy(out->Phi, in->Phi, in->N*sizeof(long_complex_t));
  out->map = in->map;
  out->time = in->time;
}

void ffluid_data_init_copy(data_ptr in, data_ptr out) {
  out->N = in->N;
  ffluid_init_data(out);
  ffluid_data_copy(in, out);
  gfluid_setup_grid(out->map, out->N);
  //ffluid_setup_grid(out);
}


void ffluid_data_fma(long_double_t op1, data_ptr op2, data_ptr op3, data_ptr out) {
  unsigned long N = op2->N;
  for (unsigned long j = 0; j < N; j++) {
    out->Q[j] = fmaq(op1, op2->Q[j], op3->Q[j]);
    out->V[j] = fmaq(op1, op2->V[j], op3->V[j]);
    //out->Q[j] = op1*op2->Q[j] + op3->Q[j];
    //out->V[j] = op1*op2->V[j] + op3->V[j];
  }
}
