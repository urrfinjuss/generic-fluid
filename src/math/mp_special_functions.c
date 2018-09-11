#include "gfluid.h"

void mpjacobi_am(mpfr_t *out, mpfr_t *in, mpfr_t *inkc) {
  /* Evaluate the Jacobi amplitude function is the easiest way to 
   * compute the elliptic sn, cn and dn simultaneously. We evaluate 
   * the Jacobi amplitude by expansion of am(u, k) at kc = 0 and 
   * evaluating the series up to kc^14. 
   *
   * The Jacobi amplitude is accurate up to 18 digits for kc = 1/20
   * 			     accurate up to 32 digits for kc = 1/40  */
  //__float128 c[25];
  //__float128 x, t;
  //__float128 se1, sh1, ch1, th1;
  //__float128 ch2;
  //__float128 ch3, sh3;
  /*
  __float128 ch4;
  __float128 sh5, ch5;
  __float128 sh7, ch7;
  */
  mpfr_t	c[25];
  mpfr_t 	x, t;
  mpfr_t 	se1, sh1, ch1, th1;
  mpfr_t	ch2;
  mpfr_t	ch3, sh3;

  mpfr_t 	mp_tmp;
  mpfr_t	mp_pi;

  mpfr_inits2(MPFR_PREC, c[0], c[1], c[2], NULL);
  mpfr_inits2(MPFR_PREC, x, t, NULL);
  mpfr_inits2(MPFR_PREC, se1, sh1, ch1, th1, NULL);
  mpfr_inits2(MPFR_PREC, ch2, ch3, sh3, NULL);
  mpfr_inits2(MPFR_PREC, mp_tmp, mp_pi, NULL);
  //mpfr_inits2(256, mp_c0, mp_c1, mp_tmp, mp_out, mp_t, NULL);
  mpfr_const_pi(mp_pi, MPFR_RNDN);
  mpfr_set(x, *in, MPFR_RNDN);
  //x = *in;
  /* quad precision trigonometry */
  /*
  se1 = 1.Q/coshq(x);
  sh1 = sinhq(x);
  ch1 = coshq(x);
  th1 = tanhq(x);
  ch2 = coshq(2.Q*x);
  ch3 = coshq(3.Q*x);
  sh3 = sinhq(3.Q*x);
  ch4 = coshq(4.Q*x);
  sh5 = sinhq(5.Q*x);
  ch5 = coshq(5.Q*x);
  sh7 = sinhq(7.Q*x);
  ch7 = coshq(7.Q*x);
  */
  /* some MPFR trigonometry */
  mpfr_set(mp_tmp, x, MPFR_RNDN);
  mpfr_sech(se1, mp_tmp, MPFR_RNDN);
  mpfr_sinh_cosh(sh1, ch1, mp_tmp, MPFR_RNDN);
  mpfr_div(th1, sh1, ch1, MPFR_RNDN);
  mpfr_mul_ui(mp_tmp, x, 2, MPFR_RNDN);
  mpfr_cosh(ch2, mp_tmp, MPFR_RNDN);
  mpfr_mul_ui(mp_tmp, x, 3, MPFR_RNDN);
  mpfr_sinh_cosh(sh3, ch3, mp_tmp, MPFR_RNDN);

  /* MPFR for the c0 coefficient */
  /* O(1) */
  mpfr_exp(mp_tmp, x, MPFR_RNDN); 
  mpfr_atan(c[0], mp_tmp, MPFR_RNDN);
  mpfr_mul_ui(c[0], c[0], 2, MPFR_RNDN);
  mpfr_div_si(mp_tmp, mp_pi, -2, MPFR_RNDN);
  mpfr_add(c[0], c[0], mp_tmp, MPFR_RNDN);
  /* MPFR for the c1 coefficient */
  /* O(xc^2) */
  mpfr_fms(c[1], x, se1, sh1, MPFR_RNDN);
  mpfr_div_ui(c[1], c[1], 4, MPFR_RNDN);
  /* O(xc^4) */
  mpfr_mul_ui(c[2], x, 2, MPFR_RNDN);  		//  2x
  mpfr_set_ui(mp_tmp, 5, MPFR_RNDN);   		//  5
  mpfr_fma(c[2], c[2], th1, mp_tmp, MPFR_RNDN); //  2x*th1 + 5
  mpfr_mul(c[2], c[2], x, MPFR_RNDN);		//  x*(5 + 2x*th1)
  mpfr_mul(c[2], c[2], se1, MPFR_RNDN);		//  x*se1*(5 + 2x*th1)
  mpfr_mul_ui(mp_tmp, sh1, 9, MPFR_RNDN);	//  9*sh1
  mpfr_sub(c[2], mp_tmp, c[2], MPFR_RNDN);	//  9*sh1 - x*se1*(5 + 2x*th1)
  mpfr_mul_si(mp_tmp, x, -4, MPFR_RNDN);	//  -4x
  mpfr_fma(c[2], mp_tmp, ch1, c[2], MPFR_RNDN); //  -4x*ch1 + 9*sh1 - x*se1*(5+2x*th1)
  mpfr_div_ui(c[2], c[2], 64, MPFR_RNDN);	//  c2/64
  /* series expansion near k = 1  */
  mpfr_pow_ui(t, *inkc, 2, MPFR_RNDN);
  /* O(xc^4) */
  /*
  c[2] = -4.Q*x*ch1 + 9.Q*sh1 - x*se1*(5.Q + 2.Q*x*th1);
  c[2] = c[2]/64.Q;
  */
  /* O(xc^6) */
  /*
  c[3] = -84.Q*x*ch1 + 8.Q*powq(x*se1, 3) + (151.Q + 12.Q*x*x - ch2)*sh1;
  c[3] += -2.Q*x*se1*(33.Q + 2.Q*x*x + 15.Q*x*th1);
  c[3] = -c[3]/1536.Q;
  */
  /* O(xc^8) */
  /*
  c[4] =  -8.Q*x*ch1*(4.Q*x*x + 285.Q) + 48.Q*powq(x*se1,3)*(x*th1 + 5.Q);
  c[4] += -x*se1*(8*x*x*x*th1 + 120.Q*x*x + 678.Q*x*th1 + 1407.Q);
  c[4] += 456.Q*x*x*sh1 + 12.Q*x*ch3 + 3747.Q*sh1 - 24.Q*sh3;
  c[4] = c[4]/49152.Q;
  */
  /* O(xc^10) */
  /*
  c[5] =     3.Q*(81595*sh1 - 740.Q*sh3 + sh5) - 60.Q*x*(64.Q*x*x + 2619)*ch1;
  c[5] += -768.Q*powq(x*se1, 5) + 160.Q*x*x*(x*x + 228.Q)*sh1;
  c[5] +=   80.Q*powq(x*se1, 3)*(8.Q*x*x + 60.Q*x*th1 + 207.Q);
  c[5] +=   -4.Q*x*x*se1*(8.Q*x*x*x + 5.Q*(40.Q*x*x + 2067.Q)*th1 + 2070.Q*x);
  c[5] += -360.Q*x*x*sh3 + 1740.Q*x*ch3 - 82740.Q*x*se1;
  c[5] = -c[5]/3932160.Q;
  */
  /* O(xc^12) */
  /*
  c[6] = -1920.Q*powq(x*se1,5)*(2.Q*x*th1 + 15.Q) + 720.Q*x*ch3*(3.Q*x*x+76.Q);
  c[6] +=   -6.Q*x*ch1*(32.Q*x*x*x*x + 17720.Q*x*x + 551175.Q);
  c[6] +=   15.Q*sh1*(464.Q*x*x*x*x + 54624.Q*x*x - 153.Q*ch2*(16.Q*x*x+49.Q) + 24.Q*ch4 + 327639.Q);
  c[6] +=   -2.Q*x*se1*(15.Q*(40.Q*x*x*x*x + 5704.Q*x*x + 51669.Q) + x*th1*(16.Q*x*x*x*x + 9780.Q*x*x + 397305.Q));
  c[6] +=  240.Q*powq(x*se1,3)*(100.Q*x*x + 1426.Q + x*th1*(8.Q*x*x + 489.Q)) - 90.Q*x*ch5;
  c[6] = c[6]/94371840.Q;
  */
  /* O(xc^14) */
  /*
  c[7] =   2520.Q*x*ch3*(5633.Q + 468.Q*x*x) - 84.Q*x*ch1*(7850715.Q + 299160.Q*x*x + 1088.Q*x*x*x*x);
  c[7] += -3780.Q*x*(15.Q*ch5 + 74668.Q*se1) + 184320.Q*powq(x*se1, 7);
  c[7] +=   112.Q*x*x*sh1*(-90.Q*ch2*(18.Q*x*x + 1193.Q) + 4.Q*(4.Q*x*x*x*x + 4350.Q*x*x + 381915.Q) + 225.Q*ch4);
  c[7] +=   224.Q*powq(x*se1,3)*(208.Q*x*x*x*x + 28200.Q*x*x + 180.Q*x*th1*(20.Q*x*x + 627.Q) + 297945.Q);
  c[7] +=   -16.Q*x*x*se1*(x*(16.Q*x*x*x*x + 19740.Q*x*x + 2085615.Q) + 105.Q*th1*(8.Q*x*x*x*x + 2508.Q*x*x + 87831.Q));
  c[7] +=-53760.Q*powq(x*se1, 5)*(4.Q*x*x + 30.Q*x*th1 + 141.Q);
  c[7] +=   -45.Q*(-21462693.Q*sh1 + 285838.Q*sh3 - 1484.Q*sh5 + sh7);
  c[7] = -c[7]/21139292160.Q;
  */
  /* O(xc^16) */
  /*
  c[8] = -322560.Q*powq(x*se1,7)*(x*th1 + 10.Q);
  c[8] +=    -84.Q*x*x*(-45.Q*sh3*(528.Q*x*x + 16051.Q) + 4.Q*sh1*(104.Q*x*x*x*x + 56095.Q*x*x + 3965880.Q) + 3750.Q*sh5);
  c[8] +=  13440.Q*powq(x*se1,5)*(280.Q*x*x + 2.Q*x*th1*(10.Q*x*x + 639.Q) + 4797.Q);
  c[8] +=    -56.Q*powq(x*se1,3)*(520.Q*x*x*(28.Q*x*x + 1845.Q) + 3.Q*x*th1*(208.Q*x*x*x*x + 51120.Q*x*x + 1134675.Q) + 8402805.Q);
  c[8] +=      1.Q*x*x*se1*(140.Q*x*(32.Q*x*x*x*x + 19188.Q*x*x + 1680561.Q) + th1*(64.Q*powq(x,6) + 143136.Q*powq(x,4)+31770900.Q*x*x + 1003691115.Q));
  c[8] =     -2.Q*c[8] - 4.Q*x*ch1*(8.Q*x*x*(32.Q*x*x*x*x + 66864.Q*x*x + 12742275.Q)+2377525815.Q);
  c[8] +=   252.Q*x*ch3*(864.Q*x*x*x*x + 121320.Q*x*x + 982075.Q);
  c[8] +=    -2.Q*1680.Q*x*ch5*(25.Q*x*x + 489.Q);
  c[8] +=   315.Q*x*(4.Q*ch7 - 11998869.Q*se1);
  c[8] +=   315.Q*(43341737.Q*sh1 - 8.Q*(81499*sh3 - 588.Q*sh5 + sh7));
  c[8] = c[8]/338228674560.Q;
  */
  /* O(xc^18) -- incomplete
   * O(xc^20) -- incomplete
   * O(xc^22) -- incomplete
   * O(xc^24) -- incomplete
   * A couple more orders would be much better, there is a
   * problem with cn(u,kc) at the ends of the interval.           */
  /* sum the terms starting with the smallest by Horner's method */
  /*
  *out = c[8];
  for (int j = 7; j > -1; j-- ) {
    // *out = -(*out)*t + c[j];
    *out = fmaq(*out, -t, c[j]);
  }
  */
  // mpfr test
  /* *out = c[1]; */
  mpfr_set(*out, c[2], MPFR_RNDN);
  for (int j = 1; j > -1; j-- ) {
    // *out = -(*out)*t + c[j];
    mpfr_fms(*out, *out, t, c[j], MPFR_RNDN);
    mpfr_neg(*out, *out, MPFR_RNDN);
  }
  mpfr_clears(c[0], c[1], c[2], NULL);
  mpfr_clears(x, t, NULL);
  mpfr_clears(se1, sh1, ch1, th1, NULL);
  mpfr_clears(ch2, ch3, sh3, NULL);
  mpfr_clears(mp_tmp, mp_pi, NULL);
}



void mp_jacobi_am2(mpfr_t *out, mpfr_t *in, mpfr_t *inkc) {
  const int			nterms = 128;
  int 				i, ii, l;
  mpfr_t			a, b, c, d, emc, u;
  mpfr_t			sn, cn, dn;
  mpfr_t			x, y;
  mpfr_t			em[nterms], en[nterms];
  mpfr_t			ca;

  /* Initialize MPFR variables */
  mpfr_inits2(MPFR_PREC, sn, cn, dn, ca, x, y, NULL);
  mpfr_inits2(MPFR_PREC, a, b, c, d, emc, u, NULL);
  for (int j = 0; j < nterms; j++) mpfr_inits2(MPFR_PREC, em[j], en[j], NULL);
  mpfr_set_float128(ca, 1.0E-27Q, MPFR_RNDN);

  /* set copies of the variables */
  mpfr_sqr(emc, *inkc, MPFR_RNDN);
  mpfr_set(u, *in, MPFR_RNDN);
  
  if ( mpfr_cmp_si(emc, 0) == 0) {
    mpfr_sech(cn, u, MPFR_RNDN);
    mpfr_set(dn, cn, MPFR_RNDN);
    mpfr_tanh(sn, u, MPFR_RNDN);
  } else {

    mpfr_set_ui(a, 1, MPFR_RNDN);   // a = 1
    mpfr_set_ui(dn, 1, MPFR_RNDN);  // dn = 1
    l = nterms-1;

    for (i = 0; i < nterms; i++) {
      mpfr_set(em[i], a, MPFR_RNDN);   	// m_array[i] = a
      mpfr_sqrt(emc, emc, MPFR_RNDN);  	// m_comp = sqrt(m_comp)
      mpfr_set(en[i], emc, MPFR_RNDN); 	// n_array[i] = m_comp
      mpfr_add(c, a, emc, MPFR_RNDN);	// c = a + m_comp
      mpfr_div_ui(c, c, 2, MPFR_RNDN);  // c = c/2
      
      mpfr_sub(x, a, emc, MPFR_RNDN);		// x = a - m_comp
      mpfr_abs(x, x, MPFR_RNDN);    		// x = |a - m_comp|
      mpfr_mul(y, ca, a, MPFR_RNDN);  		// y = CA*a
      if (mpfr_cmp(x, y) <= 0) {
        l = i;					// l = i
      	break;					// break
      }
      mpfr_mul(emc, emc, a, MPFR_RNDN);		// m_comp = m_comp * a
      mpfr_set(a, c, MPFR_RNDN);		// a = c
    }
    mpfr_mul(u, u, c, MPFR_RNDN);		// u_copy = c * u_copy
    mpfr_sin_cos(sn, cn, u, MPFR_RNDN);		// sn = sin(u_copy) 
    						// cn = cos(u_copy)
    if (mpfr_cmp_si(sn, 0) != 0 ) {		// if (sn != 0)
      mpfr_div(a, cn, sn, MPFR_RNDN);		// a = cn/sn
      mpfr_mul(c, c, a, MPFR_RNDN);		// c = c * a

      for (ii = l; 0 <= ii; ii--) {		// 
        mpfr_set(b, em[ii], MPFR_RNDN);		// b = m_array[i]
	mpfr_mul(a, a, c, MPFR_RNDN);		// a = c * a
	mpfr_mul(c, c, dn, MPFR_RNDN);		// c = c * dn

        mpfr_add(x, b, a, MPFR_RNDN);		// x = b + a
	mpfr_add(y, en[ii], a, MPFR_RNDN);	// y = n_array[i] + a
	mpfr_div(dn, y, x, MPFR_RNDN);		// dn = ( n_array[i] + a ) / ( b + a )
	mpfr_div(a, c, b, MPFR_RNDN);		// a  = c / b
      }
      mpfr_mul(a, c, c, MPFR_RNDN);		// a = c * c
      mpfr_add_ui(a, a, 1, MPFR_RNDN);		// a = c * c + 1
      mpfr_sqrt(a, a, MPFR_RNDN);		// a = sqrt( c * c + 1 )
      mpfr_ui_div(a, 1, a, MPFR_RNDN);		// a = 1 / sqrt ( c * c + 1 )

      if (mpfr_cmp_si(sn, 0) > 0) mpfr_set(sn, a, MPFR_RNDN);
      else mpfr_neg(sn, a, MPFR_RNDN);

      mpfr_mul(cn, c, sn, MPFR_RNDN); 		// cn = c * sn
    }
  } 
  //printf("sn = %39.32Qe\n", mpfr_get_float128(sn, MPFR_RNDN));
  //printf("cn = %39.32Qe\n", mpfr_get_float128(cn, MPFR_RNDN));


  mpfr_atan2(*out, sn, cn, MPFR_RNDN);
  mpfr_clears(sn, cn, dn, x, y, NULL);
  mpfr_clears(a, b, c, d, emc, u, NULL);

  for (int j = nterms-1; 0 <= j; j--) mpfr_clears(em[j], en[j], NULL);
}

void mp_jacobi_elliptick(mpfr_t *out, mpfr_t *ink) {
  mpfr_t 	one, tmp;
  
  mpfr_inits2(MPFR_PREC, one, tmp, NULL);
  mpfr_set_ui(one, 1, MPFR_RNDN);

  //mpfr_fms(tmp, *ink, *ink, one, MPFR_RNDN);
  //mpfr_neg(tmp, tmp, MPFR_RNDN);
  //mpfr_sqrt(tmp, tmp, MPFR_RNDN);

  mpfr_agm(*out, one, *ink, MPFR_RNDN);
  
  mpfr_const_pi(tmp, MPFR_RNDN);
  mpfr_div(*out, tmp, *out, MPFR_RNDN);
  mpfr_div_ui(*out, *out, 2, MPFR_RNDN);
  
  mpfr_clears(one, tmp, NULL);
}






