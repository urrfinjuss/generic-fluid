

/* declared in mapping/halfplane.c and mapping/disc.c*/
extern void ffluid_setup_grid(data_ptr in);

typedef struct conformal_map {
  unsigned long		N;	// number of gridpoints
  long_double_t		*u;	// array of gridpoints on the real line
  long_double_t		*dq;    // array of dq/du at the gridpoints
  long_double_t		kc;	// elliptic modulus complement, or \chi_c
} cmap;


/* declared in mapping/mapping.c */
void set_HTmap(cmap *map);
void free_HTmap(cmap *map);

/* Global Variables */
extern cmap this_map;
extern cmap last_map;

