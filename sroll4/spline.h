#include <unistd.h>
#include <inttypes.h>

#ifndef SPLINEDEF
#define SPLINEDEF

typedef void (*findIntervalFunction) (double, double, int, int, double,
				      int *);

enum BSPLINE_TYPE { BSPLINE_UNIFORM };

typedef struct
{
  double a;
  double b;
  int64_t k;
  int64_t g;
  int64_t Nknots;
  double *knots;
  int64_t Nvals;
  double *vals;
  int64_t Ncoefs;
  double *coefs;
  double *dp, *dm;
  int64_t type;
  int64_t verbose;
} BSpline;

BSpline *bspline_alloc (int k, int gp2);
void bspline_free (BSpline * bspline);

void bspline_init_uniform (BSpline * bspline, double a, double b);

void bspline_find_interval_uniform (double a, double b, int k, int g,
				    double x, int *j);

int bspline_value (BSpline * bspline, double x, int *istart, int *iend);

#if 0
int bspline_fit (BSpline * bspline, int Ndata, double *X, double *Y);
#endif
int bspline_eval (BSpline * bspline, int Ndata, double *X, double *Y);

int bspline_write_to_file (BSpline * spline, int filedes);

BSpline * bspline_create_from_file (int filedes);

int bspline_set_verbose(BSpline * bspline, int verbose);

int bspline_get_verbose(BSpline * bspline);

int bspline_set_default_verbose(int verbose);

int bspline_get_default_verbose(void);


#endif
