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
    fprintf(fh, "%.16QE\t", PIl*(2.0L*j/N - 1.L));
    fprintf(fh, "%.16QE\t%.16QE\n", creall(in[j]), cimagl(in[j]));
  }
  fclose(fh);
}


