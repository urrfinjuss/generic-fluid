

/* declared in mapping/halfplane.c and mapping/disc.c*/
/*
extern void ffluid_setup_grid(data_ptr in);
*/
typedef struct conformal_map {
  long_double_t		*u;	// array of gridpoints on the real line
  long_double_t		*dq;    // array of dq/du at the gridpoints
  long_double_t		kc;	// elliptic modulus complement, or \chi_c
} cmap, *cmap_ptr;


/* declared in mapping/mapping.c */
//void set_HTmap(cmap *map);
void gfluid_setup_grid(cmap *map, unsigned long N);
void free_HTmap(cmap *map);


