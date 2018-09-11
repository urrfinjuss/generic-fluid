#include "gfluid.h"

static aux_data AuxLocal;
static fft_list FFTLocal;
static __float128	qc, kc;
static __float128	Gravity;
//static __float128	Sigma;
static unsigned long 	N, kmax;

void ffluid_alloc_operators() {
  N  = DataCurr.N;
  ffluid_alloc_aux_array(&AuxLocal, 8, N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
  kc = (DataCurr.map)->kc;
  /* Construct Multipliers */
  __float128 	K0, K1;
  elliptic_k(&K0, &kc);
  elliptic_kc(&K1, &kc);
  qc = M_PIq*K1/K0;
  kmax = (unsigned long) 32*floorq(logq(10)/qc) + 1;
  Gravity = Control.Gravity;
  printf("The vertical box extends to qc = %39.32Qe\n", qc);
  printf("Positive Fourier modes are cut at k_max = %lu\n", kmax);
}

void filter_high(data_ptr in) {
  /**/
  __float128 	tmp;
  memcpy(AuxLocal.X[0], in->Q, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[1], in->V, N*sizeof(long_complex_t));
  fftwq_execute(FFTLocal.fp[0]);
  fftwq_execute(FFTLocal.fp[1]);
  for (long int j = 1; j < N/2; j++) {
    tmp = 1.Q/(1.Q + expq(powq(j - 1.44Q*kmax, 1)));
    /**/
    /* the negative modes */
    AuxLocal.Y[0][j] = AuxLocal.Y[0][j]/N;
    AuxLocal.Y[1][j] = AuxLocal.Y[1][j]/N;
    /* the positive modes */
    AuxLocal.Y[0][N-j] = AuxLocal.Y[0][N-j]*tmp/N;
    AuxLocal.Y[1][N-j] = AuxLocal.Y[1][N-j]*tmp/N;
    /**/
  }
  AuxLocal.Y[0][0] = AuxLocal.Y[0][0]/N;
  AuxLocal.Y[1][0] = AuxLocal.Y[1][0]/N;
  AuxLocal.Y[0][N/2] = AuxLocal.Y[0][N/2]/N;
  AuxLocal.Y[1][N/2] = AuxLocal.Y[1][N/2]/N;
  
  fftwq_execute(FFTLocal.bp[0]);
  fftwq_execute(FFTLocal.bp[1]);
  memcpy(in->Q, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(in->V, AuxLocal.X[1], N*sizeof(long_complex_t));
  /**/
  /* also project the nonzero modes afterward: */
  //gfluid_filter_projector(in->Q, in->Q);
  //gfluid_filter_projector(in->Q, in->Q);
  //gfluid_filter_projector(in->V, in->V);

}

void gfluid_filter_projector(long_complex_t *in, long_complex_t *out) {
  memcpy(AuxLocal.X[0], in, N*sizeof(long_complex_t));
  fftwq_execute(FFTLocal.fp[0]);
  for (long int j = 1; j < N/2; j++) {
    AuxLocal.Y[2][j]   = (cexpq(2.Q*j*qc)*AuxLocal.Y[0][j] - AuxLocal.Y[0][N-j])/2.Q/sinhq(2.Q*j*qc)/N; 
    AuxLocal.Y[2][N-j] = (AuxLocal.Y[0][j] -  cexpq(-2.Q*j*qc)*AuxLocal.Y[0][N-j])/2.Q/sinhq(2.Q*j*qc)/N;      
    /* try mpfr version of this */
  }
  AuxLocal.Y[2][0] = 0.Q;
  for (long int j = N/2; j > 0; j--) {
    AuxLocal.Y[2][0] += -powq(1.Q, j+1)*(AuxLocal.Y[0][j] - AuxLocal.Y[0][N-j])/2.Q/sinhq(1.Q*j*qc)/N;
  }
  AuxLocal.Y[2][0] += 0.5Q*AuxLocal.Y[0][0]/N;
  fftwq_execute(FFTLocal.bp[2]);
  memcpy(out, AuxLocal.X[2], N*sizeof(long_complex_t));
}

void gfluid_projector(long_complex_t *in, long_complex_t *out) {
  memcpy(AuxLocal.X[0], in, N*sizeof(long_complex_t));
  fftwq_execute(FFTLocal.fp[0]);
  for (long int j = 1; j < N/2; j++) {
    AuxLocal.Y[2][j]   = (cexpq(2.Q*j*qc)*AuxLocal.Y[0][j] - AuxLocal.Y[0][N-j])/2.Q/sinhq(2.Q*j*qc)/N; 
    AuxLocal.Y[2][N-j] = (AuxLocal.Y[0][j] -  cexpq(-2.Q*j*qc)*AuxLocal.Y[0][N-j])/2.Q/sinhq(2.Q*j*qc)/N;      
    /* try mpfr version of this */
  }
  AuxLocal.Y[2][0] = 0.Q;
  for (long int j = N/2; j > 0; j--) {
    AuxLocal.Y[2][0] += -powq(1.Q, j+1)*(AuxLocal.Y[0][j] - AuxLocal.Y[0][N-j])/2.Q/sinhq(1.Q*j*qc)/N;
  }
  AuxLocal.Y[2][0] += 0.5Q*AuxLocal.Y[0][0]/N;
  fftwq_execute(FFTLocal.bp[2]);
  memcpy(out, AuxLocal.X[2], N*sizeof(long_complex_t));
}

void gfluid_call_rhs(data_ptr in, data_ptr out) {
  cmap_ptr	map = in->map;

  filter_high(in);
  memcpy(AuxLocal.X[6], in->Q, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[7], in->V, N*sizeof(long_complex_t));
  for (long int j = 0; j < N; j++) {
    AuxLocal.X[3][j] = 2.Q*crealq(cpowq(in->Q[j],2)*conjq(in->V[j]))/(map->dq[j]);
    AuxLocal.X[4][j] = crealq((in->V[j])*conjq(in->V[j]));
  }
  gfluid_projector(AuxLocal.X[3], AuxLocal.X[3]);
  gfluid_projector(AuxLocal.X[4], AuxLocal.X[4]);
  for (long int j = 0; j < N; j++) {
    AuxLocal.X[3][j] = (map->dq[j])*AuxLocal.X[3][j];
  }
  memcpy(AuxLocal.X[5], AuxLocal.X[3], N*sizeof(long_complex_t));
  for (long int j = 0; j < N; j++) {
    AuxLocal.X[0][j] = (in->Q[j])/N;
    AuxLocal.X[1][j] = (in->V[j])/N;
    AuxLocal.X[2][j] = (AuxLocal.X[4][j])/N;
    AuxLocal.X[3][j] = (AuxLocal.X[3][j])/N;
  }
  fftwq_execute(FFTLocal.fp[0]);
  fftwq_execute(FFTLocal.fp[1]);
  fftwq_execute(FFTLocal.fp[2]);
  fftwq_execute(FFTLocal.fp[3]);
  
  AuxLocal.Y[0][0] = 0.Q;
  AuxLocal.Y[1][0] = 0.Q;
  AuxLocal.Y[2][0] = 0.Q;
  AuxLocal.Y[3][0] = 0.Q;
  for (long int j = 1; j < N/2; j++) {
    AuxLocal.Y[0][j] = -1.IQ*j*AuxLocal.Y[0][j];
    AuxLocal.Y[1][j] = -1.IQ*j*AuxLocal.Y[1][j];
    AuxLocal.Y[2][j] = -1.IQ*j*AuxLocal.Y[2][j];
    AuxLocal.Y[3][j] = -1.IQ*j*AuxLocal.Y[3][j];
    /* derivative of the positive modes too */
    AuxLocal.Y[0][N-j] = 1.IQ*j*AuxLocal.Y[0][N-j];
    AuxLocal.Y[1][N-j] = 1.IQ*j*AuxLocal.Y[1][N-j];
    AuxLocal.Y[2][N-j] = 1.IQ*j*AuxLocal.Y[2][N-j];
    AuxLocal.Y[3][N-j] = 1.IQ*j*AuxLocal.Y[3][N-j];
  }
  AuxLocal.Y[0][N/2] = 0.Q;
  AuxLocal.Y[1][N/2] = 0.Q;
  AuxLocal.Y[2][N/2] = 0.Q;
  AuxLocal.Y[3][N/2] = 0.Q;

  fftwq_execute(FFTLocal.bp[0]);
  fftwq_execute(FFTLocal.bp[1]);
  fftwq_execute(FFTLocal.bp[2]);
  fftwq_execute(FFTLocal.bp[3]);

  for (long int j = 0; j < N; j++) {
    out->Q[j] = 1.IQ*(AuxLocal.X[5][j]*AuxLocal.X[0][j] -  0.5Q*AuxLocal.X[6][j]*AuxLocal.X[3][j]);
    out->V[j] = 1.IQ*(AuxLocal.X[5][j]*AuxLocal.X[1][j] - cpowq(AuxLocal.X[7][j], 2)*AuxLocal.X[2][j]);
    out->V[j] += Gravity*(cpowq(AuxLocal.X[6][j], 2)/(map->dq[j]) - 1.Q);
  }
  filter_high(out);
  //ffluid_write_raw(out, "rhs_000.txt");
  //exit(1);
}











