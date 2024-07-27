extern int totiter;
extern double RMSD(int n, double (*motif1)[3], double (*motif2)[3]);
extern double NormDiff(int n, double (*motif1)[3], double (*motif2)[3]);
extern void motifproj(int n, double (*motif)[3], double (*in)[3], double (*out)[3]);
extern double normalize(double u[3]);
extern void mult(double a[3][3],double x[3],double ax[3]);
extern int axes(double sym[3][3],double u[3],double v[3]);
extern void cross(double u[3],double v[3],double uv[3]);
