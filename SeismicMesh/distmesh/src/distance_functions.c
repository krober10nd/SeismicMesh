// Copyright (C) 2004-2012 Per-Olof Persson. See COPYING.TXT for details.

#include <math.h>

#define szint int

// LAPACK
#ifdef __cplusplus
extern "C" {
#endif
  void dhseqr_(char*,char*,szint*,szint*,szint*,double*,szint*,double*,
               double*,void*,szint*,double*,void*,szint*);
#ifdef __cplusplus
}
#endif

// Quick routines
static inline double sqr(double x)
  { return x*x; }
static inline double length(double x, double y)
  { return sqrt(x*x+y*y); }
static inline double dot2(double const *x, double const *y)
  { return x[0]*y[0]+x[1]*y[1]; }
static inline double fmaxn(double const *x, int n) {
  double val = x[0]; int i;
  for (i=1; i<n; i++)
    if (val < x[i])
      val = x[i];
  return val;
}


// Returns non-zero if an error occurred.
static int roots(double *pol,double *rr,double *ri,szint n)
{
  double *H, *work;
  szint o=1;
  szint info;
  char chE='E',chN='N';
  H=(double*)malloc(n*n*sizeof(double));
  work=(double*)malloc(n*sizeof(double));

  memset(H,0,n*n*sizeof(double));
  int i;
  for (i=0; i<n-1; i++)
    H[1+(n+1)*i]=1.0;
  for (i=0; i<n; i++)
    H[n*i]=-pol[i+1]/pol[0];

  dhseqr_(&chE,&chN,&n,&o,&n,H,&n,rr,ri,0,&n,work,&n,&info);

  free(work); free(H);
  return info;
}

/**
 * Distance Functions
 */

double dellipse(double x0,double y0,double a,double b)
{
  double t1,t2,t4,t6,t8,t9,t11,t15,t16;
  t1=a*a; t2=b*b; t4=y0*y0; t6=x0*x0; t8=t1*t1;
  t9=t2*t1; t11=t2*t2; t15=t8*t2; t16=t11*t1;
  double pol[5]={ 1.0, 2.0*t1+2.0*t2, -t2*t4-t1*t6+t8+4.0*t9+t11,
                  -2.0*t9*t6-2.0*t9*t4+2.0*t15+2.0*t16,
                  -t16*t6+t8*t11-t15*t4 };

  double rr[4],ri[4];
  if (roots(pol,rr,ri,4) != 0)
    return NAN;

  double t=fmaxn(rr, 4);
  double x=sqr(a)*x0/(t+sqr(a));
  double y=sqr(b)*y0/(t+sqr(b));
  double d=t*sqrt(sqr(x)/sqr(sqr(a))+sqr(y)/sqr(sqr(b)));

  return d;
}

double dellipsoid(double x0,double y0,double z0,double a,double b,double c)
{
  // Begin Maple Generated
  double t1,t2,t3,t5,t6,t7,t8,t9,t10,t11,t14,t15,
         t16,t17,t18,t20,t22,t24,t41,t42,t60;

  t1=a*a; t2=b*b; t3=c*c; t5=x0*x0; t6=t1*t5; t7=t1*t1; t8=t1*t2; t9=4.0*t8;
  t10=t2*t2; t11=t1+t2; t14=t3*t3; t15=y0*y0; t16=t2*t15; t17=z0*z0;
  t18=t3*t17; t20=t7*t2; t22=t1*t10; t24=t7+t9+t10; t41=t7*t10; t42=t20+t22;
  t60=t41+4.0*t42*t3+t24*t14-t18*t7-4.0*t18*t8-t18*t10-t16*t7-
      4.0*t16*t1*t3-t16*t14-t6*t10-4.0*t6*t2*t3-t6*t14;
  double pol[7]= {
    1.0,
    2.0*t1+2.0*t2+2.0*t3,
    -t6+t7+t9+t10+4.0*t11*t3+t14-t16-t18,
    2.0*t20+2.0*t22+2.0*t24*t3+2.0*t11*t14-2.0*t16*t1-2.0*t16*t3-2.0*t6*t2-2.0*t6*t3-2.0*t18*t1-2.0*t18*t2,
    t60,
    2.0*t41*t3+2.0*t42*t14-2.0*t18*t20-2.0*t18*t22-2.0*t16*t7*t3-2.0*t16*t1*t14-2.0*t6*t10*t3-2.0*t6*t2*t14,
    t41*t14-t6*t10*t14-t18*t41-t16*t7*t14
  };
  // End Maple Generated

  double rr[6],ri[6];
  if (roots(pol,rr,ri,6) != 0)
    return NAN;

  double t=fmaxn(rr, 6);
  double x=sqr(a)*x0/(t+sqr(a));
  double y=sqr(b)*y0/(t+sqr(b));
  double z=sqr(c)*z0/(t+sqr(c));
  double d=t*sqrt(sqr(x)/sqr(sqr(a))+sqr(y)/sqr(sqr(b))+sqr(z)/sqr(sqr(c)));

  return d;
}

double dsegment(double x0,double y0,double p1x,double p1y,double p2x,double p2y)
{
  double v[2]={p2x-p1x,
               p2y-p1y};
  double w[2]={x0-p1x,
               y0-p1y};

  double c1=v[0]*w[0] + v[1]*w[1]; // dot2(v,w);
  double c2=v[0]*v[0] + v[1]*v[1]; // dot2(v,v);

  if (c1<=0)
    return length(x0-p1x,y0-p1y);
  else if (c1>=c2)
    return length(x0-p2x,y0-p2y);
  else
    return length(x0-(p1x+c1/c2*v[0]),
                  y0-(p1y+c1/c2*v[1]));
}
