#include "gfluid.h"


void elliptic_k(__float128 *out, __float128 *inkc) {
/*  Evaluate the complete elliptic integral of the 
 *  first kind by evaluating its series expansion 
 *  at k = 1 (kc = 0). 
 *
 *  Current accuracy is at least 16 digits at kc = 1/2
 *                      at least 32 digits at kc = 1/4
 *                      				*/
  __float128 c[26];
  __float128 t, xc;
  __float128 l2, lx;

  /* series expansion near k = 1  */
  t  = powq(*inkc, 2);
  xc = sqrtq(t);
  l2 = logq(2.Q);
  lx = logq(xc);
  /* O(1) */
  c[0] = 2.Q*l2 - lx;
  /* O(xc^2) */
  c[1] = -1.Q + 2.Q*l2 - lx;
  c[1] = -c[1]/4.Q;
  /* O(xc^4) */
  c[2] = -7.Q + 12.Q*l2 - 6.Q*lx;
  c[2] = 3.Q*c[2]/128.Q;
  /* O(xc^6) */
  c[3] = -37.Q + 60.Q*l2 - 30.Q*lx;
  c[3] = -5.Q*c[3]/1536.Q;
  /* O(xc^8) */
  c[4] = -533.Q + 840.Q*l2 - 420.Q*lx;
  c[4] = 35.Q*c[4]/196608.Q;
  /* O(xc^10) */
  c[5] = -1627.Q + 2520.Q*l2 - 1260.Q*lx;
  c[5] = -63.Q*c[5]/1310720.Q;
  /* O(xc^12) */
  c[6] = -18107.Q + 27720.Q*l2 - 13860.Q*lx;
  c[6] = 77.Q*c[6]/20971520.Q;
  /* O(xc^14) */
  c[7] = -237371.Q + 360360.Q*l2 - 180180.Q*lx;
  c[7] = -143.Q*c[7]/587202560.Q;
  /* O(xc^16) */
  c[8] = -95549.Q + 144144.Q*l2 - 72072.Q*lx;
  c[8] = 32175.Q*c[8]/60129542144.Q;
  /* O(xc^18) */
  c[9] = -1632341.Q + 2450448.Q*l2 - 1225224.Q*lx;
  c[9] = -60775.Q*c[9]/2164663517184.Q;
  /* O(xc^20) */
  c[10] = -155685007.Q + 232792560.Q*l2 - 116396280.Q*lx;
  c[10] = 46189.Q*c[10]/173173081374720.Q;
  /* O(xc^22) */
  c[11] = -156188887.Q + 232792560.Q*l2 - 116396280.Q*lx;
  c[11] = -29393.Q*c[11]/120946279055360.Q;
  /* O(xc^24) */
  c[12] = -3602044091.Q + 5354228880.Q*l2 - 2677114440.Q*lx;
  c[12] = 676039.Q*c[12]/69665056735887360.Q;
  /* O(xc^26) */
  c[13] = -18051406831.Q + 26771144400.Q*l2 - 13385572200.Q*lx;
  c[13] = -1300075.Q*c[13]/724516590053228544.Q;
  /* O(xc^28) */
  c[14] = -7751493599.Q + 11473347600.Q*l2 - 5736673800.Q*lx;
  c[14] = 5014575.Q*c[14]/1288029493427961856.Q;
  /* O(xc^30) */
  c[15] = -225175759291.Q + 332727080400.Q*l2 - 166363540200.Q*lx;
  c[15] = -646323.Q*c[15]/5152117973711847424.Q;
  /* O(xc^32) */
  c[16] = -13981692518567.Q + 20629078984800.Q*l2 - 10314539492400.Q*lx;
  c[16] = 20036013.Q*c[16]/10551537610161863524352.Q;
  /* O(xc^34) */
  c[17] = -14000078506967.Q + 20629078984800.Q*l2 - 10314539492400.Q*lx;
  c[17] = -116680311.Q*c[17]/65227687044636974514176.Q;
  /* O(xc^36) */
  c[18] = -98115155543129.Q + 144403552893600.Q*l2 - 72201776446800.Q*lx;
  c[18] = 756261275*c[18]/3130928978142574776680448.Q;
  /* O(xc^38) */
  c[19] = -3634060848592973.Q + 5342931457063200.Q*l2 - 2671465728531600.Q*lx;
  c[19] = -1472719325.Q*c[19]/237950602338835683027714048.Q;
  /* O(xc^40) */
  c[20] = -3637485804655193.Q + 5342931457063200.Q*l2 - 2671465728531600.Q*lx;
  c[20] =  2297442147.Q*c[20]/390483039735525223430094848.Q;
  /* O(xc^42) */
  c[21] = -149264130644602513.Q + 219060189739591200.Q*l2 - 109530094869795600.Q*lx;
  c[21] = -4485482287.Q*c[21]/32800575337784118768127967232.Q;
  /* O(xc^44) */
  c[22] = -6423336258393807859.Q + 9419588158802421600.Q*l2 - 4709794079401210800.Q*lx;
  c[22] =  17534158031.Q*c[22]/5772901259450004903190522232832.Q;
  /* O(xc^46) */
  c[23] = -6427886784074388739.Q + 9419588158802421600.Q*l2 - 4709794079401210800.Q*lx;
  c[23] = -514589420475.Q*c[23]/177035638623133483697842681806848.Q;
  /* O(xc^48) */
  c[24] = -302306920271471321183.Q + 442720643463713815200.Q*l2 - 221360321731856907600.Q*lx;
  c[24] =  2687300306925.Q*c[24]/45321123487522171826647726542553088.Q;
  /* sum the terms starting with the smallest by Horner's method */
  *out = c[24];
  for (int j = 23; j > -1; j-- ) {
    //*out = -(*out)*t + c[j];
    *out = fmaq(*out, -t, c[j]);
  }
}

void jacobi_am(__float128 *out, __float128 *in, __float128 *inkc) {
  /* Evaluate the Jacobi amplitude function is the easiest way to 
   * compute the elliptic sn, cn and dn simultaneously. We evaluate 
   * the Jacobi amplitude by expansion of am(u, k) at kc = 0 and 
   * evaluating the series up to kc^14. 
   *
   * The Jacobi amplitude is accurate up to 18 digits for kc = 1/20
   * 			     accurate up to 32 digits for kc = 1/40  */
  __float128 c[25];
  __float128 x, t;
  __float128 se1, sh1, ch1, th1;
  __float128 ch2;
  __float128 ch3, sh3;
  __float128 ch4;
  __float128 sh5, ch5;
  __float128 sh7, ch7;

  x = *in;
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
  /* series expansion near k = 1  */
  t = powq(*inkc, 2);
  /* O(1) */
  c[0] = -0.5Q*PIq + 2.Q*atanq(expq(x));
  /* O(xc^2) */
  c[1] = -x*se1 + sh1;
  c[1] = -0.25Q*c[1];
  /* O(xc^4) */
  c[2] = -4.Q*x*ch1 + 9.Q*sh1 - x*se1*(5.Q + 2.Q*x*th1);
  c[2] = c[2]/64.Q;
  /* O(xc^6) */
  c[3] = -84.Q*x*ch1 + 8.Q*powq(x*se1, 3) + (151.Q + 12.Q*x*x - ch2)*sh1;
  c[3] += -2.Q*x*se1*(33.Q + 2.Q*x*x + 15.Q*x*th1);
  c[3] = -c[3]/1536.Q;
  /* O(xc^8) */
  c[4] =  -8.Q*x*ch1*(4.Q*x*x + 285.Q) + 48.Q*powq(x*se1,3)*(x*th1 + 5.Q);
  c[4] += -x*se1*(8*x*x*x*th1 + 120.Q*x*x + 678.Q*x*th1 + 1407.Q);
  c[4] += 456.Q*x*x*sh1 + 12.Q*x*ch3 + 3747.Q*sh1 - 24.Q*sh3;
  c[4] = c[4]/49152.Q;
  /* O(xc^10) */
  c[5] =     3.Q*(81595*sh1 - 740.Q*sh3 + sh5) - 60.Q*x*(64.Q*x*x + 2619)*ch1;
  c[5] += -768.Q*powq(x*se1, 5) + 160.Q*x*x*(x*x + 228.Q)*sh1;
  c[5] +=   80.Q*powq(x*se1, 3)*(8.Q*x*x + 60.Q*x*th1 + 207.Q);
  c[5] +=   -4.Q*x*x*se1*(8.Q*x*x*x + 5.Q*(40.Q*x*x + 2067.Q)*th1 + 2070.Q*x);
  c[5] += -360.Q*x*x*sh3 + 1740.Q*x*ch3 - 82740.Q*x*se1;
  c[5] = -c[5]/3932160.Q;
  /* O(xc^12) */
  c[6] = -1920.Q*powq(x*se1,5)*(2.Q*x*th1 + 15.Q) + 720.Q*x*ch3*(3.Q*x*x+76.Q);
  c[6] +=   -6.Q*x*ch1*(32.Q*x*x*x*x + 17720.Q*x*x + 551175.Q);
  c[6] +=   15.Q*sh1*(464.Q*x*x*x*x + 54624.Q*x*x - 153.Q*ch2*(16.Q*x*x+49.Q) + 24.Q*ch4 + 327639.Q);
  c[6] +=   -2.Q*x*se1*(15.Q*(40.Q*x*x*x*x + 5704.Q*x*x + 51669.Q) + x*th1*(16.Q*x*x*x*x + 9780.Q*x*x + 397305.Q));
  c[6] +=  240.Q*powq(x*se1,3)*(100.Q*x*x + 1426.Q + x*th1*(8.Q*x*x + 489.Q)) - 90.Q*x*ch5;
  c[6] = c[6]/94371840.Q;
  /* O(xc^14) */
  c[7] =   2520.Q*x*ch3*(5633.Q + 468.Q*x*x) - 84.Q*x*ch1*(7850715.Q + 299160.Q*x*x + 1088.Q*x*x*x*x);
  c[7] += -3780.Q*x*(15.Q*ch5 + 74668.Q*se1) + 184320.Q*powq(x*se1, 7);
  c[7] +=   112.Q*x*x*sh1*(-90.Q*ch2*(18.Q*x*x + 1193.Q) + 4.Q*(4.Q*x*x*x*x + 4350.Q*x*x + 381915.Q) + 225.Q*ch4);
  c[7] +=   224.Q*powq(x*se1,3)*(208.Q*x*x*x*x + 28200.Q*x*x + 180.Q*x*th1*(20.Q*x*x + 627.Q) + 297945.Q);
  c[7] +=   -16.Q*x*x*se1*(x*(16.Q*x*x*x*x + 19740.Q*x*x + 2085615.Q) + 105.Q*th1*(8.Q*x*x*x*x + 2508.Q*x*x + 87831.Q));
  c[7] +=-53760.Q*powq(x*se1, 5)*(4.Q*x*x + 30.Q*x*th1 + 141.Q);
  c[7] +=   -45.Q*(-21462693.Q*sh1 + 285838.Q*sh3 - 1484.Q*sh5 + sh7);
  c[7] = -c[7]/21139292160.Q;
  /* O(xc^16) */
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
  /* O(xc^18) -- incomplete
   * O(xc^20) -- incomplete
   * O(xc^22) -- incomplete
   * O(xc^24) -- incomplete
   * A couple more orders would be much better, there is a
   * problem with cn(u,kc) at the ends of the interval.           */
  /* sum the terms starting with the smallest by Horner's method */
  *out = c[8];
  for (int j = 7; j > -1; j-- ) {
    //*out = -(*out)*t + c[j];
    *out = fmaq(*out, -t, c[j]);
  }
}























void jacobi_sn(__float128 *out, __float128 *in, __float128 *inkc) {
  __float128 c[7];
  __float128 t, tmp;
  __float128 sec1, tan1;
  __float128 sin2, cos2;
  __float128 sin4, cos4;
  __float128 sin6;
  
  tan1 = tanhq(*in);
  sec1 = 1.Q/coshq(*in);
  sin2 = sinhq(2.Q*(*in));
  cos2 = coshq(2.Q*(*in));
  sin4 = sinhq(4.Q*(*in));
  cos4 = coshq(4.Q*(*in));
  sin6 = sinhq(6.Q*(*in));

  /* series expansion near k = 1  */
  t = powq(*inkc, 2);
  /*   O(m-1)^6   */
  c[6] =    1920.Q*(70.Q*powq(*in,2) + 6.Q*(*in)*sin2 + 1599.Q)*(*in)*cos2;
  c[6] +=     -4.Q*(*in)*(2880.Q*(4.Q*(*in)*tan1 + 29.Q)*powq((*in)*sec1,4)*sec1*sec1);
  c[6] +=    960.Q*powq((*in)*sec1,3)*sec1*(348.Q*powq(*in,2) + 2.Q*(16.Q*powq(*in,2) + 675.Q)*(*in)*tan1 + 3675.Q);
  c[6] +=     -4.Q*(*in)*powq(sec1,2)*(96.Q*powq(*in,2)*(116.Q*powq(*in,2)+6125.Q) + 2.Q*(*in)*tan1*(256.Q*powq(*in,4)+54000.Q*powq(*in,2)+930105.Q)+3161385.Q);
  c[6] += -29520.Q*(*in)*cos4 - 1920.Q*powq(*in,2)*sin2*(4.Q*powq(*in,2)+489.Q);
  c[6] +=    -45.Q*(87141.Q*sin2-892.Q*sin4 + sin6 - 384164.Q*tan1);
  c[6] = c[6]/188743680.Q;
  /*   O(m-1)^5   */
  c[5] =   -960.Q*powq((*in)*sec1, 5)*sec1;
  c[5] +=    10.Q*(*in)*cos2*(32.Q*powq(*in,2) + 1449.Q);
  c[5] +=   -45.Q*((80.Q*powq(*in,2) + 457.Q)*sin2 - 3.Q*sin4 - 2225.Q*tan1);
  c[5] +=    -1.Q*(*in)*(128.Q*powq(*in,4) + 12440.Q*powq(*in,2) + 30.Q*(64.Q*powq(*in,2) + 1415.Q)*(*in)*tan1 + 73965.Q)*powq(sec1, 2);
  c[5] +=    60.Q*powq((*in)*sec1,3)*(16.Q*powq(*in,2) + 96.Q*(*in)*tan1 + 311.Q)*sec1;
  c[5] +=   -60.Q*(*in)*cos4;
  c[5] =  -c[5]/983040.Q;
  /*   O(m-1)^4   */
  c[4] =     48.Q*powq((*in)*sec1, 3)*(4.Q*(*in)*tan1 + 19.Q)*sec1;
  c[4] +=     3.Q*((-32.Q*powq(*in, 2) - 339.Q)*sin2 + sin4 + 1899.Q*tan1);
  c[4] +=  -608.Q*powq((*in)*sec1,2)*(*in);
  c[4] +=    -4.Q*powq((*in)*sec1,2)*(16.Q*powq(*in,2) + 591.Q)*tan1;
  c[4] += -4275.Q*(*in)*sec1*sec1;
  c[4] +=   600.Q*(*in)*cos2;
  c[4] = c[4]/49152.Q;
  /*   O(m-1)^3   */
  c[3] =   12.Q*powq((*in)*sec1, 3)*sec1; 
  c[3] +=  -8.Q*powq((*in)*sec1, 2)*(*in);
  c[3] += -42.Q*powq((*in)*sec1, 2)*tan1;
  c[3] += -15.Q*sin2 +  6.Q*(*in)*cos2;
  c[3] += 105.Q*tan1 - 81.Q*(*in)*sec1*sec1;
  c[3] = -c[3]/768.Q;
  /*   O(m-1)^2   */
  c[2] = -4.Q*powq((*in)*sec1,2)*tan1 - sin2 + 11.Q*tan1 - 9.Q*(*in)*powq(sec1,2);
  c[2] = c[2]/64.Q;
  /*   O(m-1)     */
  c[1] = -0.25Q*(tan1 - (*in)*powq(sec1,2));
  /*   O(1)       */  
  c[0] = tan1;
 
  *out = 1.Q*c[6];
  *out = -(*out)*t + c[5];
  *out = -(*out)*t + c[4];
  *out = -(*out)*t + c[3];
  *out = -(*out)*t + c[2];
  *out = -(*out)*t + c[1];
  *out = -(*out)*t + c[0];
  
}

void jacobi_elliptic(__float128 *sn, __float128 *cn, __float128 *dn, __float128 *in, __float128 *inkc) {
  __float128 phi;
  jacobi_am(&phi, in, inkc);
  *sn = sinq(phi);
  *cn = cosq(phi);
  *dn = sqrtq( (*cn)*(*cn) + (*inkc)*(*inkc)*(*sn)*(*sn));
}
