

/* declared in mapping/halfplane.c and mapping/disc.c*/
/*
extern void ffluid_setup_grid(data_ptr in);
*/
typedef struct conformal_map {
  __float128		*u;	// array of gridpoints on the real line
  __float128		*dq;    // array of dq/du at the gridpoints
  __float128		kc;	// elliptic modulus complement, or \chi_c
} cmap, *cmap_ptr;


/* declared in mapping/mapping.c */
//void set_HTmap(cmap *map);
void gfluid_setup_grid(cmap *map, unsigned long N);
void free_HTmap(cmap *map);


