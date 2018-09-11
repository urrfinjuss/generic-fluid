/* declared in memory/memory.c */
extern void ffluid_memory_module();
extern void ffluid_init_data(data_ptr in);
extern void ffluid_alloc_aux_array(aux_data_ptr in, unsigned long NArrays, unsigned long NElements);
extern void ffluid_alloc_fft_plans(aux_data_ptr in, fft_list_ptr out);
extern void ffluid_dealloc_aux_array(aux_data_ptr in);


/* declared in memory/array_func.c  */
extern void gfluid_data_init_copy(data_ptr in, data_ptr out);
extern void gfluid_data_copy(data_ptr in, data_ptr out);
extern void gfluid_data_fma(long_double_t op1, data_ptr op2, data_ptr op3, data_ptr out);
