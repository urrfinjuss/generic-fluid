#include "gfluid.h"
#include <math.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_sf_ellint.h>


void evaluate_ht_transform(long double XiC, long N) {
  double 	sn2, cn2, dn2;
  double 	sn1, cn1, dn1;
  double 	s1,   c1,  d1;
  long double 	q, w, mult;
  long double 	m = sqrtl(1.0L - XiC)*sqrtl(1.0L + XiC);
  long double 	K = gsl_sf_ellint_Kcomp (m, 2e-32L)/PIl;


  printf("K = %.16Le\n", K);
  FILE *fh = fopen("ht_gsl_elliptic.txt", "w");
  fprintf(fh, "# 1. q 2. u\n# XiC = %.12Le\n\n", XiC);
  //fprintf(fh, "%23.16Le\t%23.16Le\n", -PIl, -PIl);
  q = 2.0L*K*PIl/N;
  q = PIl*K;
  // GSL
  gsl_sf_elljac_e(q, m*m, &s1, &c1, &d1);
  // Mathematica
  //s1 = 0.0144350857170492880983762481717183761288228605124106;
  //c1 = 0.9998958087222595412253205358662337872247913655774692;
  //d1 = 0.9998968506887788099188653100774316915829981478902475;

  sn1 = s1; cn1 = c1; dn1 = d1;
  
  printf("k = %.16Le\n", m);
  printf("m = k^2 = %.16Le\n", m*m);
  printf("q = %.16Le\n", q);
  printf("sn(q, k) = %.17e\n", sn1);
  printf("cn(q, k) = %.17e\n", cn1);
  printf("dn(q, k) = %.17e\n", dn1);
  printf("w(q) = %.17Le\n", 2.0L*asinl(XiC*sn1/dn1));
  fprintf(fh, "%23.16Le\t%23.16Le\n", 0.L, 0.L);
  for (long j = 0; j < N/2; j++) {
    //q = PIl*(2.0L*j/N - 1.0L);
    //gsl_sf_elljac_e(2.*K*q, m*m, &sn2, &cn2, &dn2);
    //gsl_sf_elljac_e(K*q, m*m, &sn1, &cn1, &dn1);
    //mult = (dn2 + 1.)/(XiC*XiC + dn2 + cn2);
    //mult = 1./(1 - m*m*sn1*sn1);
    w = 2.0L*asinl(XiC*sn1/dn1);
    fprintf(fh, "%23.16Le\t%23.16Le\n", PIl*(2.0L*(j+1)/N), w);
    //naive double angle
    cn2 =     (cn1*c1 - sn1*s1*dn1*d1)/(1.L - m*m*sn1*sn1*s1*s1);
    sn2 =     (sn1*c1*d1 + s1*cn1*dn1)/(1.L - m*m*sn1*sn1*s1*s1);
    dn2 = (dn1*d1 - m*m*sn1*s1*cn1*c1)/(1.L - m*m*sn1*sn1*s1*s1);
  
    sn1 = sn2;
    cn1 = cn2;
    dn1 = dn2;
  }
  fclose(fh);
}



void evaluate_jacobi_sn() { 
  __float128 u_loc, l_loc, m_loc, k_loc, sn;
  u_loc = 0.5Q*PIq;

  /*  Xc = 1e-1  */
  l_loc = 1.0E-5Q;
  m_loc = 1.0Q - l_loc*l_loc;
  k_loc = sqrtq(m_loc);
 
  char sn_true[] = "9.188134412891826335550473828773884328468797723720536795487983050e-1";

  double sn_gsl, cn_gsl, dn_gsl;
  double ud = 0.5*PIl;
  double ld = 1.0E-5;
  double md = 1.0 - ld*ld;

  gsl_sf_elljac_e(ud, md, &sn_gsl, &cn_gsl, &dn_gsl);
  printf("u   =\t%39.32Qe\n", u_loc);
  printf("k =\t%39.32Qe\n", k_loc);
  jacobi_sn(&sn, &u_loc, &l_loc);
  printf("GSL sn:\tsn(u,k) = %23.17e\n", sn_gsl);
  printf("Exact:\tsn(u,k) = %s\n", sn_true);
  printf("Series:\tsn(u,k) = %38.32Qe\t(sn expansion)\n", sn);
  
  //char ek_true[] = "3.695637362989874677809954195262550984801e+0";
  char ek_true[] = "2.156515647499643235438674998800322028864e+0";
  __float128 ek_series;
  elliptic_k(&ek_series, &l_loc);
  printf("Exact:\tK(Xc = 1/2) = %s\n", ek_true);
  printf("Series:\tK(Xc = 1/2) = %38.32Qe\n", ek_series);
  
  char am_true[] = "1.1870500303181522121928284444791e+0";
  __float128 am_series;

  jacobi_am(&am_series, &u_loc, &l_loc);
  printf("Exact:\tam(u,k) = %s\n", am_true);
  printf("Series:\tam(u,k) = %38.32Qe\n", am_series);

  __float128 cn, dn;
  jacobi_elliptic(&sn, &cn, &dn, &u_loc, &l_loc);
  printf("Series:\tsn(u,k) = %38.32Qe\t(am expansion)\n", sn);
  printf("GSL:   \tsn(u,k) = %38.32e\t\n", sn_gsl);
  printf("Series:\tcn(u,k) = %38.32Qe\t(am expansion)\n", cn);
  printf("GSL:   \tcn(u,k) = %38.32e\t\n", cn_gsl);
  printf("Series:\tdn(u,k) = %38.32Qe\t(am expansion)\n", dn); 
  printf("GSL:   \tdn(u,k) = %38.32e\t\n", dn_gsl);


  
  printf("\n\n");

  /*  Xc = 1e-2  */
  l_loc = 1.0E-2Q;
  m_loc = 1.0Q - l_loc*l_loc;
  k_loc = sqrtq(m_loc);

  char sn_true2[] = "9.171690263699496985112238317155966318713510569951891319754256808e-1";

  ld = 1.0E-2;
  md = 1.0 - ld*ld;
  gsl_sf_elljac_e(ud, md, &sn_gsl, &cn_gsl, &dn_gsl);
  printf("u   =\t%39.32Qe\n", u_loc);
  printf("k =\t%39.32Qe\n", k_loc);
  jacobi_sn(&sn, &u_loc, &l_loc);
  printf("Series:\tsn(u,k) = %38.32Qe\n", sn);
  printf("GSL sn:\tsn(u,k) = %23.17e\n", sn_gsl);
  printf("Exact:\tsn(u,k) = %s\n", sn_true2);
}

void elliptic_demo() {
  __float128 u, kc;
  __float128 am, sn, cn, dn;
  __float128 k0;
  
  kc = 1e-2Q;
  // find the period of elliptic functions:
  elliptic_k(&k0, &kc);
  k0 = k0/PIq; 
  u = PIq*k0; // this might be too extreme
  jacobi_am(&am, &u, &kc); 
  jacobi_elliptic(&sn, &cn, &dn, &u, &kc);
  printf("Full Report:\n");
  printf("Elliptic K(kc) = %39.32Qe\n", PIq*k0);
  printf("u = %39.32Qe\tkc = %39.32Qe\n", u, kc);
  printf("Using series expansion:\n");
  printf("am(u, k) = %39.32Qe\n", am);
  printf("sn(u, k) = %39.32Qe\n", sn);
  printf("cn(u, k) = %39.32Qe\n", cn);
  printf("dn(u, k) = %39.32Qe\n", dn);
  
  long int	N = 256;
  __float128	w, q;
  FILE *fh = fopen("ht_transform.txt","w");
  fprintf(fh, "# 1. q 2. w\n");
  fprintf(fh, "# transform type is Hale-Tee N = %ld\tkc = %39.32Qe\n\n", N, kc);
  for (long int j = 0; j < N; j++) {
    q = PIq*k0*(2.Q*j/N - 1.Q);
    jacobi_elliptic(&sn, &cn, &dn, &q, &kc);
    w = 2.Q*asinq(kc*sn/dn);
    fprintf(fh, "%39.32Qe\t%39.32Qe\n", q/k0, w);
  }
  fclose(fh);
  

}






