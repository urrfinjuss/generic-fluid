#include "gfluid.h"

/* note: remove math.h after all has been done in quadmath */
#include <math.h>


void ffluid_write_array(long_complex_t *in, unsigned long N, char *fname) {
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. u 2. X 3. Y.\n\n");
  for (unsigned long j = 0; j < N; j++) {
    fprintf(fh, "%.16LE\t", PIl*(2.0L*j/N - 1.L));
    fprintf(fh, "%.16LE\t%.16LE\n", creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}


void gfluid_write_map(char *fname) {
  char  full_path[160];
  __float128	q;

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. q 2. u(q) 3. q_u\n");
  fprintf(fh, "# Hale-Tee:\tkc = %.4Qe\tN = %lu\n\n", this_map.kc, this_map.N);
  for (long int j = 0; j < this_map.N; j++) {
    q = PIq*(2.Q*j/this_map.N - 1.Q);
    fprintf(fh, "%39.32Qe\t", q);
    fprintf(fh, "%39.32Qe\t", this_map.u[j]);
    fprintf(fh, "%39.32Qe\n", this_map.dq[j]);
  }
  fclose(fh);

}
