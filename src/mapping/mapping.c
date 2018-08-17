#include "gfluid.h"

void set_HTmap(cmap *map) {
  unsigned long N = map->N;
  __float128 sn, cn, dn;
  __float128 k0, kc, q;

  kc = map->kc;
  map->u  = fftwq_malloc(N*sizeof(__float128));
  map->dq = fftwq_malloc(N*sizeof(__float128));
  /*  set the HT map  */
  elliptic_k(&k0, &kc);  
  k0 = k0/PIq;
  for (long int j = 0; j < N; j++) {
    q = PIq*k0*(2.Q*j/N - 1.0Q);
    jacobi_elliptic(&sn, &cn, &dn, &q, &kc);
    map->u[j] = 2.Q*asinq(kc*sn/dn);
    map->dq[j] = 0.5Q*dn/(k0*kc);
  }
}

void free_HTmap(cmap *map) {
  fftwq_free(map->u);
  fftwq_free(map->dq);
}

