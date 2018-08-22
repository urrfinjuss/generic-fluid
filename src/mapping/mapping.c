#include "gfluid.h"

void gfluid_setup_grid(cmap_ptr map, unsigned long N) {
  __float128 sn, cn, dn, am;
  __float128 k0, q;

  map->kc = Control.transform_kc;  
  elliptic_k(&k0, &(map->kc));  
  jacobi_elliptic(&sn, &cn, &dn, &q, &(map->kc));
  k0 = k0/PIq;
  map->u[0] = -PIq;
  map->dq[0] = 0.5Q/k0;
  for (long int j = 1; j < N/2; j++) {
    //q = PIq*k0*(2.Q*j/N - 1.0Q);
    q = 2.Q*PIq*k0*j/N;
    jacobi_elliptic(&sn, &cn, &dn, &q, &(map->kc));
    jacobi_am(&am, &q, &(map->kc));
    // the gridpoints
    map->u[j]   = -PIq + 2.Q*am;
    map->u[N-j] = -1.0Q*map->u[j];
    // the jacobians
    map->dq[j]   = 0.5Q/(k0*dn);
    map->dq[N-j] = map->dq[j];
  }
  map->u[N/2] = 0;
  map->dq[N/2] = 0.5Q/(map->kc)/k0;
}

