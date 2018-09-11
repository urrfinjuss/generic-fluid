#include "gfluid.h"

void gfluid_setup_grid(cmap_ptr map, unsigned long N) {
  __float128 	k0; 
  mpfr_t 	x, y, mp_l, mp_am, mp_theta, mp_kc;
  mpfr_t	mp_u, mp_q, mp_tmp, mp_pi;
  mpfr_t	K;

  map->kc = Control.transform_kc;  
  //elliptic_k(&k0, &(map->kc)); 
  mpfr_inits2(MPFR_PREC, mp_kc, K, NULL);
  mpfr_set_float128(mp_kc, map->kc, MPFR_RNDN);
  mp_jacobi_elliptick(&K, &mp_kc);
  k0 = mpfr_get_float128(K, MPFR_RNDN);
  //printf("By series (quad):\t%39.32Qe\n", k0);
  //printf("By agm (MPFR):\t\t%39.32Qe\n", mpfr_get_float128(K, MPFR_RNDN));
  mpfr_clears(K, NULL);

  //jacobi_elliptic(&sn, &cn, &dn, &q, &(map->kc));
  k0 = k0*M_2_PIq;
  map->u[0] = -M_PIq;
  map->dq[0] = 1.0Q/k0;
 
  mpfr_inits2(MPFR_PREC, x, y, mp_l, mp_am, mp_theta, NULL);
  mpfr_inits2(MPFR_PREC, mp_q, mp_tmp, mp_u, mp_pi, K, NULL);
  //mpfr_set_float128(mp_kc, map->kc, MPFR_RNDN);
  mpfr_set_float128(mp_tmp, M_PIq*k0, MPFR_RNDN);
  mpfr_const_pi(mp_pi, MPFR_RNDN);

  for (long int j = 1; j < N/2; j++) {
    /* MPFR q */
    mpfr_mul_si(mp_q, mp_tmp, j, MPFR_RNDN);
    mpfr_div_ui(mp_q, mp_q, N, MPFR_RNDN);
    //mpjacobi_am(&mp_am, &mp_q, &mp_kc);  // by series
    mp_jacobi_am2(&mp_am, &mp_q, &mp_kc);  // by recipes
    /*
    if (j == N/2 - 1) {
        //mpfr_set_ui(mp_kc, 0, MPFR_RNDN);
        mpjacobi_am(&mp_am, &mp_q, &mp_kc);
    	printf("By series:\t%39.32Qe\n", mpfr_get_float128(mp_am, MPFR_RNDN));
	mp_jacobi_am2(&mp_am, &mp_q, &mp_kc);
    	printf("By recipes:\t%39.32Qe\n", mpfr_get_float128(mp_am, MPFR_RNDN));
	exit(1);
    }
    */
    // the gridpoints (MPFR)
    mpfr_mul_ui(mp_u, mp_am, 2, MPFR_RNDN);
    mpfr_sub(mp_u, mp_u, mp_pi, MPFR_RNDN);
    map->u[j] = mpfr_get_float128(mp_u, MPFR_RNDN);
    map->u[N-j] = -(map->u[j]);
    // geometric form (MPFR)		VII
    mpfr_sin_cos(x, y, mp_am, MPFR_RNDN);
    mpfr_div(y, y, mp_kc, MPFR_RNDN);
    mpfr_atan2(mp_theta, x, y, MPFR_RNDN);
    mpfr_sin(y, mp_theta, MPFR_RNDN);
    mpfr_div(y, y, mp_kc, MPFR_RNDN);
    mpfr_div(mp_l, x, y, MPFR_RNDN);
    // set in quad precision 
    map->dq[j] = mpfr_get_float128(mp_l, MPFR_RNDN);
    map->dq[j] = 1.0Q/k0/(map->dq[j]);
    map->dq[N-j] = map->dq[j];
  }
  map->u[N/2] = 0.Q;
  map->dq[N/2] = 1.0Q/(map->kc)/k0;

  mpfr_clears(x, y, mp_theta, mp_l, mp_am, mp_kc, mp_q, mp_tmp, K, NULL);
}

