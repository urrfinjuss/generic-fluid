#include "gfluid.h"

static sim_data SimLocal;
static aux_data AuxLocal;
static fft_list FFTLocal;
//static long_complex_t z[5], beta[5];

void ffluid_alloc_equations() {
  gfluid_data_init_copy(&DataCurr, &SimLocal);
  ffluid_alloc_aux_array(&AuxLocal, 5, DataCurr.N);
  ffluid_alloc_fft_plans(&AuxLocal, &FFTLocal);
}

void gfluid_make_spectrum(data_ptr in) {
  unsigned long N = in->N;
  memcpy(AuxLocal.X[0], in->Q, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[1], in->V, N*sizeof(long_complex_t));
  fftwq_execute(FFTLocal.fp[0]);
  fftwq_execute(FFTLocal.fp[1]);
  memcpy(in->Q, AuxLocal.Y[0], N*sizeof(long_complex_t));
  memcpy(in->V, AuxLocal.Y[1], N*sizeof(long_complex_t));
}

void gfluid_get_natural_variables(data_ptr in) {
  /* find the natural variables: the map and the potential */
  gfluid_get_surface(in);
  //get_potential(in);

  /* find the constants of motion */
  //find_energy(in);
  //find_momentum(in);
}

void gfluid_get_surface(data_ptr in) {
  unsigned long 	N = in->N;
  //cmap_ptr		map = in->map;
  __float128		y0 = 0.Q;
  //long_complex_t	y1 = 0.Q;

  memcpy(AuxLocal.X[0], in->Q, N*sizeof(long_complex_t));
  for (long int j = 0; j < N; j++) {
    AuxLocal.X[0][j] = cimagq(cpowq(in->Q[j],-2))/N;
    AuxLocal.X[1][j] = crealq(cpowq(in->Q[j],-2))/N;
    AuxLocal.X[2][j] =       (cpowq(in->Q[j],-2) - 1.Q)/N;  // this has z 2 tildas!
    //AuxLocal.X[2][j] =       (cpowq(in->Q[j],-2) - 1.Q/(map->dq[j]))/N;  // this has z 1 tildas!
  }
  fftwq_execute(FFTLocal.fp[0]);
  fftwq_execute(FFTLocal.fp[1]);
  fftwq_execute(FFTLocal.fp[2]);
  for (long int j = 1; j < N/2; j++) {
    /* get antiderivative of Zq */
    AuxLocal.Y[2][j]   =  1.0IQ*AuxLocal.Y[2][j]/j;
    AuxLocal.Y[2][N-j] = -1.0IQ*AuxLocal.Y[2][N-j]/j;
    /* get antiderivative of Yq */
    AuxLocal.Y[0][j]   =  1.0IQ*AuxLocal.Y[0][j]/j;
    AuxLocal.Y[0][N-j] = -1.0IQ*AuxLocal.Y[0][N-j]/j;
  }
  AuxLocal.Y[0][N/2] = 0.Q;
  AuxLocal.Y[0][0]   = 0.Q;
  for (long int j = N/2-1; j > 0; j--) {
    y0 += crealq( AuxLocal.Y[0][  j]*conjq(AuxLocal.Y[1][  j]));
    y0 += crealq( AuxLocal.Y[0][N-j]*conjq(AuxLocal.Y[1][N-j]));
  }
  AuxLocal.Y[2][N/2] = 0.Q;
  AuxLocal.Y[2][0] = -1.0IQ*crealq(y0);
  fftwq_execute(FFTLocal.bp[2]);
  memcpy(in->Z, AuxLocal.X[2], N*sizeof(long_complex_t));
  /* Mean level computation debugging */
  /* */
  long_double_t		ml = 0.Q;
  //printf("Zero mode y0 is %39.32Qe\n", y0);
  for (long int j = 0; j < N; j++) {
    ml += cimagq(in->Z[j])*crealq(cpowq(in->Q[j], -2));
  }
  ml = 2.Q*PIq*ml/N;
  //printf("Mean Level at %39.32Qe\n", ml);
  /* */
}

/*
void ffluid_call_rhs(data_ptr in, data_ptr out) {
  unsigned long 	N = in->N;
  long_complex_t	w1 = cexpl( 1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	w2 = cexpl(-1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	b1U = 0.0L, b2U = 0.0L;

  if (in->l == 1.0L) {
    w2 = 0.0L;
    w1 = 0.0L;
  }
  memcpy(AuxLocal.X[0], in->R, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[1], in->V, N*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.fp[0]); 
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);

  for (unsigned long j = 0; j < N/2; j++) {
    AuxLocal.Y[0][j] = -1.IL*j*AuxLocal.Y[0][j]/N;
    AuxLocal.Y[1][j] = -1.IL*j*AuxLocal.Y[1][j]/N;
  }
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j]= 2.0L*creall(in->V[j]*in->du[j]*conjl(in->R[j]*in->R[j]))/N;
    AuxLocal.X[3][j]= (in->V[j]*conjl(in->V[j]) + 2.0L*Control.Sigma*(in->R[j]*conjl(in->R[j]) + 2.0L*cimagl(AuxLocal.X[0][j]*conjl(in->R[j]))))/N;
  }
  fftwl_execute(FFTLocal.fp[2]);
  fftwl_execute(FFTLocal.fp[3]);
  AuxLocal.Y[4][0] = 0.0L;
  b1U = AuxLocal.Y[2][N/2+1];
  b2U = AuxLocal.Y[2][N/2-1];
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    b1U = b1U*w1 + AuxLocal.Y[2][N/2+1+j];
    b2U = b2U*w2 + AuxLocal.Y[2][N/2-1-j];
  }
  memset(AuxLocal.Y[2]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[3]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[4]+N/2, 0, N/2*sizeof(long_complex_t));
  AuxLocal.Y[2][0] = 0.5L*AuxLocal.Y[2][0];
  */
  /* compute the proper zero mode of the equation: */
  /* z_t = i U z_u                                 */
  /* to plot the surface correctly                 */
  /*
  // derivative via memmove
  memmove(AuxLocal.Y[3], AuxLocal.Y[3]+1, (N-1)*sizeof(long_complex_t));
  for (unsigned long j = 0; j < N/2; j++) {
    // this is derivative wrt to u
    // AuxLocal.Y[3][j] = -1.0IL*j*AuxLocal.Y[3][j]; 
    // this is derivative wrt to \xi (verified)
    AuxLocal.Y[3][j] = -1.0L*(j+1)*AuxLocal.Y[3][j];
    AuxLocal.Y[4][j] = -1.0IL*j*AuxLocal.Y[2][j];
  }
  fftwl_execute(FFTLocal.bp[2]);
  fftwl_execute(FFTLocal.bp[3]);
  fftwl_execute(FFTLocal.bp[4]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j] += 0.5L*(b1U*w1 - b2U*w2);
  }
  for (long int j = 0; j < N; j++) { 
    */
    /* these are the evolution equations for S,V */
    //AuxLocal.X[0][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[0][j] - AuxLocal.X[4][j]*in->R[j] + 1.0IL*AuxLocal.X[2][j]*in->R[j])/N;
    //AuxLocal.X[1][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j] - AuxLocal.X[3][j]*in->R[j])/N;
    /* this is the evolution equation for Q,V */
    /*
    AuxLocal.X[0][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[0][j] - 0.5L*AuxLocal.X[4][j]*in->R[j] + 0.5IL*AuxLocal.X[2][j]*in->R[j])/N;
    AuxLocal.X[1][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j] - AuxLocal.X[3][j]*in->R[j]*in->R[j])/N;
    //AuxLocal.X[0][j] *= in->du[j];
    //AuxLocal.X[1][j] *= in->du[j];
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[1]+N/2, 0, N/2*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
  */
    /*
    out->Q[j] = 0.5IL*(2.L*AuxLocal.X[0][j]*AuxLocal.X[2][j]-AuxLocal.X[4][j]*in->Q[j]);
    out->V[j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j]-in->Q[j]*in->Q[j]*AuxLocal.X[3][j]);
    out->Q[j] = out->Q[j]*in->du[j];
    out->V[j] = out->V[j]*in->du[j];
    */
  /*  
  ffluid_write_array(out->R, N, "rhsR.txt");
  ffluid_write_array(out->V, N, "rhsV.txt");
  exit(0);
  */
//}



// experimental
/*
void ffluid_call_rhsX(data_ptr in, data_ptr out) {
  unsigned long 	N = in->N;
  long_complex_t	w1 = cexpl( 1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	w2 = cexpl(-1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	b1U = 0.0L, b2U = 0.0L;

  if (in->l == 1.0L) {
    w2 = 0.0L;
    w1 = 0.0L;
  }
  memcpy(AuxLocal.X[0], in->R, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[1], in->V, N*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.fp[0]); 
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);

  for (unsigned long j = 0; j < N/2; j++) {
    AuxLocal.Y[0][j] = -1.IL*j*AuxLocal.Y[0][j]/N;
    AuxLocal.Y[1][j] = -1.IL*j*AuxLocal.Y[1][j]/N;
  }
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j]= 2.0L*creall(in->V[j]*in->du[j]*conjl(in->R[j]*in->R[j]))/N;
    AuxLocal.X[3][j]= (in->V[j]*conjl(in->V[j]) + 2.0L*Control.Sigma*(in->R[j]*conjl(in->R[j]) + 2.0L*cimagl(AuxLocal.X[0][j]*conjl(in->R[j]))))/N;
  }
  fftwl_execute(FFTLocal.fp[2]);
  fftwl_execute(FFTLocal.fp[3]);
  AuxLocal.Y[4][0] = 0.0L;
  b1U = AuxLocal.Y[2][N/2+1];
  b2U = AuxLocal.Y[2][N/2-1];
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    b1U = b1U*w1 + AuxLocal.Y[2][N/2+1+j];
    b2U = b2U*w2 + AuxLocal.Y[2][N/2-1-j];
  }
  memset(AuxLocal.Y[2]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[3]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[4]+N/2, 0, N/2*sizeof(long_complex_t));
  AuxLocal.Y[2][0] = 0.5L*AuxLocal.Y[2][0];
  // derivative via memmove
  // memmove(AuxLocal.Y[3], AuxLocal.Y[3]+1, (N-1)*sizeof(long_complex_t));
  for (unsigned long j = 0; j < N/2; j++) {
    // this is derivative wrt to u
    // AuxLocal.Y[3][j] = -1.0IL*j*AuxLocal.Y[3][j]; 
    // this is derivative wrt to \xi (unverified)
    AuxLocal.Y[3][j] = -1.0L*(j+1)*AuxLocal.Y[3][j+1];
    AuxLocal.Y[4][j] = -1.0IL*j*AuxLocal.Y[2][j];
  }
  fftwl_execute(FFTLocal.bp[2]);
  fftwl_execute(FFTLocal.bp[3]);
  fftwl_execute(FFTLocal.bp[4]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j] += 0.5L*(b1U*w1 - b2U*w2);
  }
  for (long int j = 0; j < N; j++) {
    //AuxLocal.X[0][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[0][j] - AuxLocal.X[4][j]*in->R[j] + 1.0IL*AuxLocal.X[2][j]*in->R[j])/N;
    //AuxLocal.X[1][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j] - AuxLocal.X[3][j]*in->R[j])/N;
    AuxLocal.X[0][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[0][j] - 0.5L*AuxLocal.X[4][j]*in->R[j] + 0.5IL*AuxLocal.X[2][j]*in->R[j])/N;
    AuxLocal.X[1][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j] - AuxLocal.X[3][j]*in->R[j]*in->R[j])/N;
    //AuxLocal.X[0][j] *= in->du[j];
    //AuxLocal.X[1][j] *= in->du[j];
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[1]+N/2, 0, N/2*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
}
*/











// experimental by text

/*
void ffluid_call_rhsXX(data_ptr in, data_ptr out) {
  unsigned long 	N = in->N;
  long_complex_t	w1 = cexpl( 1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	w2 = cexpl(-1.IL*in->u0 - 2.0L*atanhl(in->l));
  long_complex_t	b1U = 0.0L, b2U = 0.0L;

  if (in->l == 1.0L) {
    w2 = 0.0L;
    w1 = 0.0L;
  }
  memcpy(AuxLocal.X[0], in->R, N*sizeof(long_complex_t));
  memcpy(AuxLocal.X[1], in->V, N*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.fp[0]); 
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2);
  memset(AuxLocal.Y[1]+N/2, 0, N/2);

  for (unsigned long j = 0; j < N/2; j++) {
    AuxLocal.Y[0][j] = -1.IL*j*AuxLocal.Y[0][j]/N;
    AuxLocal.Y[1][j] = -1.IL*j*AuxLocal.Y[1][j]/N;
  }
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j]= 2.0L*creall(in->V[j]*conjl(in->R[j]*in->R[j]))/N;
    AuxLocal.X[3][j]= (in->V[j]*conjl(in->V[j]) + 0.5L*Control.Sigma*(in->R[j]*conjl(in->R[j]) + 2.0L*cimagl(AuxLocal.X[0][j]*conjl(in->R[j]))))/N;
  }
  fftwl_execute(FFTLocal.fp[2]);
  fftwl_execute(FFTLocal.fp[3]);
  AuxLocal.Y[4][0] = 0.0L;
  b1U = AuxLocal.Y[2][N/2+1];
  b2U = AuxLocal.Y[2][N/2-1];
  for (unsigned long j = 1; j < N/2 - 1; j++) {
    b1U = b1U*w1 + AuxLocal.Y[2][N/2+1+j];
    b2U = b2U*w2 + AuxLocal.Y[2][N/2-1-j];
  }
  memset(AuxLocal.Y[2]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[3]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[4]+N/2, 0, N/2*sizeof(long_complex_t));
  AuxLocal.Y[2][0] = 0.5L*AuxLocal.Y[2][0];

  for (unsigned long j = 0; j < N/2; j++) {
    AuxLocal.Y[3][j] = -1.0IL*j*AuxLocal.Y[3][j]; 
    AuxLocal.Y[4][j] = -1.0IL*j*AuxLocal.Y[2][j];
  }
  fftwl_execute(FFTLocal.bp[2]);
  fftwl_execute(FFTLocal.bp[3]);
  fftwl_execute(FFTLocal.bp[4]);
  for (unsigned long j = 0; j < N; j++) {
    AuxLocal.X[2][j] += 0.5L*(b1U*w1 - b2U*w2);
  }
  for (long int j = 0; j < N; j++) {
    AuxLocal.X[0][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[0][j] - 0.5L*AuxLocal.X[4][j]*in->R[j] + 0.5IL*AuxLocal.X[2][j]*in->R[j])/N;
    AuxLocal.X[1][j] = 1.0IL*(AuxLocal.X[2][j]*AuxLocal.X[1][j] - AuxLocal.X[3][j]*in->R[j]*in->R[j] + 1.IL*AuxLocal.X[2][j]*in->V[j])/N;
  }
  fftwl_execute(FFTLocal.fp[0]);
  fftwl_execute(FFTLocal.fp[1]);
  memset(AuxLocal.Y[0]+N/2, 0, N/2*sizeof(long_complex_t));
  memset(AuxLocal.Y[1]+N/2, 0, N/2*sizeof(long_complex_t));
  fftwl_execute(FFTLocal.bp[0]);
  fftwl_execute(FFTLocal.bp[1]);
  memcpy(out->R, AuxLocal.X[0], N*sizeof(long_complex_t));
  memcpy(out->V, AuxLocal.X[1], N*sizeof(long_complex_t));
}
*/











