#include "gfluid.h"


void ffluid_write_natural(data_ptr in, char *fname) {
  unsigned long N = in->N;
  __float128	time = in->time;
  __float128 	q;
  cmap_ptr	map = in->map;
  char full_path[160];

  // 
  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. q 2. u 3. Z 5. Phi\n");
  fprintf(fh, "# Time = %.12Qe\tkc = %.16Qe\n\n", DataCurr.time, map->kc);
  for (unsigned long j = 0; j < N; j++) {
    q = PIq*(2.Q*j/N - 1);
    fprintf(fh, "%23.16QE\t%23.16QE\t", q, map->u[j]);
    fprintf(fh, "%23.16QE\t%23.16QE\t", crealq(in->Z[j]), cimagq(in->Z[j]));
    fprintf(fh, "%23.16QE\t%23.16QE\n", crealq(in->Phi[j]), cimagq(in->Phi[j]));
  }
  fclose(fh);
}

void ffluid_write_raw(data_ptr in, char *fname) {
  unsigned long N = in->N;
  __float128	time = in->time;
  __float128 	q;
  cmap_ptr	map = in->map;
  char full_path[160];

  // 
  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. q 2. u 3. Q 4. V\n");
  fprintf(fh, "# Time = %.12Qe\tkc = %.16Qe\n\n", DataCurr.time, map->kc);
  for (unsigned long j = 0; j < N; j++) {
    q = PIq*(2.Q*j/N - 1);
    fprintf(fh, "%.16QE\t%.16QE\t", q, map->u[j]);
    fprintf(fh, "%.16QE\t%.16QE\t", crealq(in->Q[j]), cimagq(in->Q[j]));
    fprintf(fh, "%.16QE\t%.16QE\n", crealq(in->V[j]), cimagq(in->V[j]));
  }
  fclose(fh);
}

void ffluid_write_spectrum(data_ptr in, char *fname) {
  unsigned long N = in->N;
  cmap_ptr	map = in->map;
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. k 2. |Q_k| 4. |V_k|\n");
  fprintf(fh, "# Time = %.12Qe\tkc = %.16Qe\n\n", DataCurr.time, map->kc);
  for (unsigned long j = 0; j < N/2; j++) {
    fprintf(fh, "%.16QE\t", -1.0L*j);
    fprintf(fh, "%.16QE\t%.16QE\n", cabsq(in->Q[j]), cabsq(in->V[j]));
  }
  fclose(fh);
}

void ffluid_start_log(char *fname) {
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. Time 2.-3. Z(-pi) 4. Volume (const) 5. KEnergy 6. SEnergy\n\n");
  fclose(fh);
}

void ffluid_append_to_log(data_ptr in, char *fname) {
  unsigned long N = in->N;
  char full_path[160];

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "a");
  fprintf(fh, "%.12Qe\t", DataCurr.time);
  fprintf(fh, "%.16Qe\t", crealq(in->Z[0]));
  fprintf(fh, "%.16Qe\t", cimagq(in->Z[0]));
  fprintf(fh, "%.16Qe\t", DataCurr.Volume);
  fprintf(fh, "%.16Qe\t", DataCurr.Hamiltonian);
  fprintf(fh, "%.16Qe\n", DataCurr.SurfaceEnergy);
  fclose(fh);
}

void gfluid_write_map(cmap_ptr in, unsigned long N, char *fname) {
  char  full_path[160];
  __float128	q;

  strcpy(full_path, Control.data_path);
  strcat(full_path, fname);
  FILE *fh = fopen(full_path, "w");
  fprintf(fh, "# 1. q 2. u(q) 3. q_u\n");
  fprintf(fh, "# Hale-Tee:\tkc = %.4Qe\tN = %lu\n\n", in->kc, N);
  for (long int j = 0; j < N; j++) {
    q = PIq*(2.Q*j/N - 1.Q);
    fprintf(fh, "%39.32Qe\t", q);
    fprintf(fh, "%39.32Qe\t", in->u[j]);
    fprintf(fh, "%39.32Qe", in->dq[j]);
    fprintf(fh, "\n");
  }
  fclose(fh);
}
