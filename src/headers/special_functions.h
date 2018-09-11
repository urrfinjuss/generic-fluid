

/* Special functions  */
extern void jacobi_am(__float128 *out, __float128 *in, __float128 *inkc);
extern void jacobi_sn(__float128 *out, __float128 *in, __float128 *inkc);
extern void jacobi_elliptic(__float128 *sn, __float128 *cn, __float128 *dn, __float128 *in, __float128 *inkc);

extern void elliptic_k(__float128 *out, __float128 *inkc);
extern void elliptic_kc(__float128 *out, __float128 *inkc);

/* MP Special functions  */
extern void mpjacobi_am(mpfr_t *out, mpfr_t *in, mpfr_t *inkc);  // by series, deprecated
extern void mp_jacobi_am2(mpfr_t *out, mpfr_t *in, mpfr_t *inkc);
extern void mp_jacobi_elliptick(mpfr_t *out, mpfr_t *ink);

