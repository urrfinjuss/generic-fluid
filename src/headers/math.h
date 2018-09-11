
/* declared in math/equations.c */
extern void gfluid_get_natural_variables(data_ptr in);
extern void gfluid_get_surface(data_ptr in);
extern void gfluid_make_spectrum(data_ptr in);

/* declared in math/equations.c */
extern void ffluid_alloc_equations();

/* declared in math/operators.c */
extern void ffluid_alloc_operators();
extern void gfluid_projector(long_complex_t *in, long_complex_t *out);
extern void gfluid_filter_projector(long_complex_t *in, long_complex_t *out);
extern void filter_high(data_ptr in);
extern void gfluid_call_rhs(data_ptr in, data_ptr out);



