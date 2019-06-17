/* Prototypes for all user accessible RANDLIB routines */


extern float genexp(float av);
extern void genmn(float *parm,float *x,float *work);
extern float gennor(float av,float sd);
extern float genunf(float low,float high);
extern void gscgn(long getset,long *g);
extern long ignlgi(void);
extern void initgn(long isdtyp);
extern long mltmod(long a,long s,long m);
extern float ranf(void);
extern void setall(long iseed1,long iseed2);
extern void setgmn(float *meanv,float *covm,long p,float *parm);
extern float sexpo(void);
extern float snorm(void);







