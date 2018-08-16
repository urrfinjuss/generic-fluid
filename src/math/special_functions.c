#include "gfluid.h"

void jacobi_sn(__float128 *out, __float128 *in, __float128 *ink) {
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
  t = (1.Q + (*ink))*(1.Q - (*ink));
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
