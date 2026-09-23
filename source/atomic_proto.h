/* atomicdata.c */
int get_atomic_data(char masterfile[]);
/* atomicdata_sub.c */
int atomicdata2file(void);
int index_lines(void);
int index_phot_top(void);
int index_inner_cross(void);
void indexx(int n, float arrin[], int indx[]);
int limit_lines(double freqmin, double freqmax);
int check_xsections(void);
int hr1985_setup_table(double logbeta[], double loggamma[], double log_g[], int npts);
double q21(struct lines *line_ptr, double t);
double q12(struct lines *line_ptr, double t);
double a21(struct lines *line_ptr);
double upsilon(int n_coll, double u0);
void skiplines(FILE *fptr, int nskip);
/* atomicdata_init.c */
int init_atomic_data(void);
