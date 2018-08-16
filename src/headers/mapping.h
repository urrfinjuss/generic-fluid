
/* declared in mapping/mapping.c */
extern void evaluate_ht_transform(long double XiC, long N);

/* declared in mapping/disc.c */
extern void ffluid_disk();
/* declared in mapping/halfplane.c */
extern void ffluid_halfplane();

/* declared in mapping/halfplane.c and mapping/disc.c*/
extern void ffluid_setup_grid(data_ptr in);
