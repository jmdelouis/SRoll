/********************************************
 *
 *   spline.c
 *
 *
 *
 ********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
// For mac os x
//#include <Accelerate/Accelerate.h>
// For others
//#include "mkl_lapack.h"

#include "spline.h"

#define MARKER_LEN 17
#define MARKER_START "**BSPLINE_START**"
#define MARKER_KNOTS "**BSPLINE_KNOTS**"
#define MARKER_VALS  "**BSPLINE_VALS***"
#define MARKER_COEFS "**BSPLINE_COEFS**"
#define MARKER_END   "**BSPLINE_END****"

static int bspline_default_verbose = 0;

/********************************************************
 *   bplsine_alloc                                      *
 *                                                      *
 *        - k : order of the spline (cubic = 3)         *
 *        - gp2 : number of knots (gp2-1 intervals)     *
 *                                                      *
 *      return: pointer to allocated bspline struct     *
 *                                                      */

BSpline *
bspline_alloc (int k, int gp2)
{
  BSpline *self = NULL;

  self = (BSpline *) calloc (1, sizeof (BSpline));
  if (!self)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "BSpline alloc failed 1\n");
      goto fail;
    }
  self->k = k;
  self->g = gp2 - 2;		/* number of interior knots */
  self->Nknots = self->g + 2 + 2 * self->k;	/* total number of knots */
  self->knots = (double *) malloc (self->Nknots * sizeof (double));
  if (!self->knots)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "knots alloc failed! (%" PRId64 ")\n", self->Nknots);
      goto fail;
    }
  self->Nvals = self->g + self->k + 1;
  self->vals = (double *) malloc (self->Nvals * sizeof (double));
  if (!self->vals)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "vals alloc failed! (%" PRId64 ")\n", self->Nvals);
      goto fail;
    }

  self->Ncoefs = self->Nvals;
  self->coefs = (double *) malloc (self->Ncoefs * sizeof (double));
  if (!self->coefs)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "coefs alloc failed! (%" PRId64 ")\n", self->Ncoefs);
      goto fail;
    }

  self->dp = (double *) malloc ((self->k + 1) * sizeof (double));
  if (!self->dp)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "dp alloc failed! (%" PRId64 ")\n", self->k + 1);
      goto fail;
    }

  self->dm = (double *) malloc ((self->k + 1) * sizeof (double));
  if (!self->dm)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "dm alloc failed! (%" PRId64 ")\n", self->k + 1);
      goto fail;
    }

  self->verbose = bspline_default_verbose;

  return self;
fail:
  bspline_free (self);
  return NULL;
}

/********************************************************
 *   bplsine_free                                       *
 *                                                      *
 *      free allocated data of a bspline                *
 *                                                      */

void
bspline_free (BSpline * bspline)
{
  int verbose;
  if (!bspline)
    {
      if (bspline_default_verbose)
        fprintf(stderr, "Nothing to deallocate\n");
      return;
    }
  verbose = bspline->verbose;
  free (bspline->knots);
  free (bspline->vals);
  free (bspline->coefs);
  free (bspline->dp);
  free (bspline->dm);
  free (bspline);
  if (verbose) 
    fprintf(stderr, "bspline deallocated\n");
}

/********************************************************
 *   bplsine_init_uniform                               *
 *                                                      *
 *      initialize an allocated bspline struct with     *
 *      uniform knots between a and b                   *
 *                                                      *
 *      - bspline : a pointer to an allocated bspline   *
 *      - a, b : interval of the bspline                *
 *                                                      */

void
bspline_init_uniform (BSpline * bspline, double a, double b)
{
  int i;

  int k, g, Nknots;
  double *knots;

  if (!bspline)
    return;

  bspline->type = BSPLINE_UNIFORM;
  bspline->a = a;
  bspline->b = b;

  k = bspline->k;
  g = bspline->g;
  knots = bspline->knots;
  Nknots = bspline->Nknots;

  for (i = 0; i < k; i++)
    knots[i] = a;

  for (i = k; i < k + g + 1; i++)
    knots[i] = a + (i - k) / (g + 1.) * (b - a);

  for (i = k + g + 1; i < Nknots; i++)
    knots[i] = b;
}

/********************************************************
 *   bplsine_find_interval_uniform                      *
 *                                                      *
 *      internal function used to fond the interval     *
 *      in which falls a given number                   *
 *                                                      *
 *      - a, b : interval of the spline                 *
 *      - k, g : parameters of the spline               *
 *      - x : value for which to find the interval      *
 *                                                      *
 *      - *j : interval number (output)                 *
 *                                                      */

void
bspline_find_interval_uniform (double a, double b, int k, int g, double x,
			       int *j)
{
  *j = -1;

  if (x < a || x > b)
    return;

  *j = (int) floor ((x - a) / (b - a) * (g + 1.)) + k;
  if (*j > g + k)
    *j = g + k;
  //*j = (*j < g + k) ? *j : g + k;
}

/********************************************************
 *   bplsine_value                                      *
 *                                                      *
 *     evaluate the b-splines at position x             *
 *     the results is put inside internal structure     *
 *                                                      *
 *      - bspline : a pointer to an initialized         *
 *                  b-spline structure                  *
 *      - x : value for which to compute the b-splines  *
 *                                                      *
 *      - *istart : first non-zero in vals (output)     *
 *      - *iend : last+1 non-zero in vals (output)      *
 *                                                      *
 *     return: 0 if problem, 1 if ok                    *
 *           the result is put inside bspline->vals     *
 *           array                                      *
 *                                                      */

int
bspline_value (BSpline * bspline, double x, int *istart, int *iend)
{
  int k, g; // Nknots, Nvals;
  double a, b;
  double *knots;
  findIntervalFunction find_interval;
  double *N, *t;
  double *dp, *dm;

  int i;

  if (!bspline)
    return 0;

  switch (bspline->type)
    {
    case BSPLINE_UNIFORM: 
      find_interval = bspline_find_interval_uniform; 
      break;
    default: 
      return 0;
    }

  k = bspline->k;
  g = bspline->g;
  a = bspline->a;
  b = bspline->b;
  knots = bspline->knots;
  //Nknots = bspline->Nknots;
  N = bspline->vals;
  //Nvals = bspline->Nvals;

  t = knots;			/* to follow de Boor's notations */

  dp = bspline->dp;
  dm = bspline->dm;

  /* Find interval */
  find_interval (a, b, k, g, x, &i);
  if (i < 0)
    {
      *istart = -1;
      *iend = -2;
      return 0;
    }

  /* interval of vals affected */
  *istart = i - k;
  *iend = i;


  N[0 + *istart] = 1.;
  for (int s = 1; s <= k; s++)
    {
      double old, prev;
      dp[s - 1] = t[i + s] - x;
      dm[s - 1] = x - t[i + 1 - s];
      prev = 0.0;
      old = N[0 + *istart];
      for (int r = 1; r <= s; r++)
	{
	  double M;
	  M = old / (dp[r - 1] + dm[s + 1 - r - 1]);
	  old = N[r + *istart];
	  N[r - 1 + *istart] = prev + dp[r - 1] * M;
	  N[r + *istart] = dm[s + 1 - r - 1] * M;
	  prev = N[r + *istart];
	}
    }

  return 1;
}

#if 0

/********************************************************
 *   bplsine_fit                                        *
 *                                                      *
 *    fit b-spline coefficients using given data points *
 *                                                      *
 *      - bspline : a pointer to an initialized         *
 *                  b-spline structure                  *
 *      - Ndata : number of data points used for fit    *
 *      - X[] : array of x data points (size=Ndata)     *
 *      - Y[] : array of y data points (size=Ndata)     *
 *                                                      *
 *     return: 0 if problem, 1 if ok                    *
 *           the result is put inside bspline->coefs    *
 *           array                                      *
 *                                                      */

int
bspline_fit (BSpline * bspline, int Ndata, double *X, double *Y)
{
  double *AB;			// Matrix AtA, (ncoefs x ncoefs), band diagonal with K-1 superdiag
  clock_t start, end;
  int Ncoefs;
  double *coefs;
  int Nvals;
  double *vals;
  int k;

  k = bspline->k;
  Ncoefs = bspline->Ncoefs;
  coefs = bspline->coefs;
  Nvals = bspline->Nvals;
  vals = bspline->vals;

  if (bspline->verbose)
    fprintf (stderr, "Entering spline_fit...\n");

  AB = (double *) calloc (Ncoefs * (k + 1), sizeof (double));
  memset (coefs, 0, Ncoefs * sizeof (double));

  if (bspline->verbose)
    fprintf (stderr, "Filling matrix and vector...\n");

  start = clock ();
  for (int r = 0; r < Ndata; r++)
    {
      int istart, iend;
      bspline_value (bspline, X[r], &istart, &iend);
      for (int i = istart; i <= iend; i++)
	{
	  double vi;
	  vi = vals[i];
	  coefs[i] += Y[r] * vi;
	  for (int i2 = i; i2 <= iend; i2++)
	    {
	      double vj;
	      vj = vals[i2];
	      AB[(k + i - i2) + i2 * (k + 1)] += vi * vj;
	    }
	}
    }
  end = clock ();
  if (bspline->verbose)
    fprintf (stderr, "Filling matrix and vector in %lf seconds\n",
             ((double) (end - start)) / CLOCKS_PER_SEC);

  // Solve the system
  char *uplo = "U";
  int ldab = k + 1;
  int ldb = Ncoefs;
  int nrhs = 1;
  int info;
  int NN = Ncoefs;

  start = clock ();
  dpbsv_ (uplo, &NN, &k, &nrhs, AB, &ldab, coefs, &ldb, &info);
  end = clock ();
  if (bspline->verbose)
    fprintf (stderr, "Solving system in %lf seconds [info=%d]\n",
             ((double) (end - start)) / CLOCKS_PER_SEC, info);

  // Free allocated memory
  free (AB);

  return info;

}
#endif

/********************************************************
 *   bplsine_eval                                       *
 *                                                      *
 *    eval spline using fit coefs at positions X        *
 *    put results in given array                        *
 *                                                      *
 *      - bspline : a pointer to an initialized         *
 *                  b-spline structure                  *
 *      - Ndata : number of points to evaluate          *
 *      - X[] : array of x data points (size=Ndata)     *
 *                                                      *
 *      - Y[] : (output) array of y data points         *
 *                       (size=Ndata)                   *
 *                                                      *
 *     return: 0 if problem, 1 if ok                    *
 *                                                      */

int
bspline_eval (BSpline * bspline, int Ndata, double *X, double *Y)
{
  int i, ii;
  int istart, iend;

  double *vals;
  double *coefs;

  vals = bspline->vals;
  coefs = bspline->coefs;

  for (i = 0; i < Ndata; i++)
    {
      bspline_value (bspline, X[i], &istart, &iend);
      Y[i] = 0.0;
      for (ii = istart; ii <= iend; ii++)
	Y[i] += vals[ii] * coefs[ii];
    }

  return 1;
}


/********************************************************
 *   bplsine_write_to_file                              *
 *                                                      *
 *    Write a BSpline to a file                         *
 *                                                      *
 *      - bspline : the bspline to write                *
 *      - f : a pointer to a FILE structure             *
 *                                                      *
 *     return: 0 if problem, 1 if ok                    *
 *                                                      */

int 
bspline_write_to_file (BSpline * bspline, int filedes)
{
  //size_t n;

  if (! bspline)
    return 0;

  write (filedes, (void *) MARKER_START, sizeof (char) * MARKER_LEN);
  write (filedes, (void *) bspline, sizeof (BSpline));

  write (filedes, (void *) MARKER_KNOTS, sizeof (char) * MARKER_LEN);
  write (filedes, (void *) bspline->knots, sizeof (double) * bspline->Nknots);

  write (filedes, (void *) MARKER_VALS, sizeof (char) * MARKER_LEN);
  write (filedes, (void *) bspline->vals, sizeof (double) * bspline->Nvals);

  write (filedes, (void *) MARKER_COEFS, sizeof (char) * MARKER_LEN);
  write (filedes, (void *) bspline->coefs, sizeof (double) * bspline->Ncoefs);

  write (filedes, (void *) MARKER_END, sizeof (char) * MARKER_LEN);

  return 1;
}


/********************************************************
 *   bplsine_create_from_file                           *
 *                                                      *
 *    Create a BSpline using data in FILE               *
 *                                                      *
 *      - f : a pointer to a FILE structure             *
 *                                                      *
 *     return: a pointer to the allocated BSpline       *
 *             or NULL if problem                       *
 *                                                      */

BSpline *
bspline_create_from_file (int filedes)
{
  //size_t n;
  BSpline *self = NULL;
  char marker[MARKER_LEN + 1];

  marker[MARKER_LEN] = '\0';

  self = (BSpline *) calloc (1, sizeof (BSpline));

  if (!self)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "BSpline alloc failed 1\n");
      goto fail;
    }

  read (filedes, (void *) marker, MARKER_LEN);
  if (strncmp (marker, MARKER_START, MARKER_LEN) != 0 )
    {
      if (bspline_default_verbose)
        fprintf (stderr, "bad file format\n");
      goto fail;          
    }
  
  read (filedes, (void *) self, sizeof (BSpline));

  self->coefs = NULL;
  self->vals = NULL;
  self->knots = NULL;
  self->dp = NULL;
  self->dm = NULL;

  /* Checking and reading KNOTS */
  read (filedes, (void *) marker, MARKER_LEN);
  if (strncmp (marker, MARKER_KNOTS, MARKER_LEN) != 0 )
    {
      if (bspline_default_verbose)
        fprintf (stderr, "bad file format [KNOTS]\n");
      goto fail;          
    }  

  self->knots = (double *) malloc (self->Nknots * sizeof (double));
  if (!self->knots)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "knots alloc failed! (%" PRId64 ")\n", self->Nknots);
      goto fail;
    }
  read (filedes, (void *) self->knots, self->Nknots * sizeof (double));

  /* Checking and reading VALS */
  read (filedes, (void *) marker, MARKER_LEN);
  if (strncmp (marker, MARKER_VALS, MARKER_LEN) != 0 )
    {
      if (bspline_default_verbose)
        fprintf (stderr, "bad file format [%s] [VALS]\n", marker);
      goto fail;          
    }  
  
  self->vals = (double *) malloc (self->Nvals * sizeof (double));
  if (!self->vals)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "vals alloc failed! (%" PRId64 ")\n", self->Nvals);
      goto fail;
    }
  read (filedes, (void *) self->vals, self->Nvals * sizeof (double));

  /* Checking and reading COEFS */
  read (filedes, (void *) marker, MARKER_LEN);
  if (strncmp (marker, MARKER_COEFS, MARKER_LEN) != 0 )
    {
      if (bspline_default_verbose)
        fprintf (stderr, "bad file format [COEFS]\n");
      goto fail;          
    }  
  
  self->coefs = (double *) malloc (self->Ncoefs * sizeof (double));
  if (!self->coefs)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "coefs alloc failed! (%" PRId64 ")\n", self->Ncoefs);
      goto fail;
    }
  read (filedes, (void *) self->coefs, self->Ncoefs * sizeof (double));

  read (filedes, (void *) marker, MARKER_LEN);
  if (strncmp (marker, MARKER_END, MARKER_LEN) != 0 )
    {
      if (bspline_default_verbose)
        fprintf (stderr, "bad file format [END]\n");
      goto fail;
    }
  
  self->dp = (double *) malloc ((self->k + 1) * sizeof (double));
  if (!self->dp)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "dp alloc failed! (%" PRId64 ")\n", self->k + 1);
      goto fail;
    }

  self->dm = (double *) malloc ((self->k + 1) * sizeof (double));
  if (!self->dm)
    {
      if (bspline_default_verbose) 
        fprintf (stderr, "dm alloc failed! (%" PRId64 ")\n", self->k + 1);
      goto fail;
    }

  self->verbose = bspline_default_verbose;
  
  return self;

 fail:
  bspline_free(self);
  return NULL;
}

int 
bspline_set_verbose (BSpline * bspline, int verbose)
{
  int old_verbose;

  if (bspline)
    {
      old_verbose = bspline->verbose;
      bspline->verbose = verbose;
      return old_verbose;
    }
  return -1;
}

int 
bspline_get_verbose (BSpline * bspline)
{
  if (bspline)
    return bspline->verbose;  
  return -1;
}

int 
bspline_set_default_verbose (int verbose)
{
  int old_verbose;

  old_verbose = bspline_default_verbose;
  bspline_default_verbose = verbose;

  return old_verbose;
}

int 
bspline_get_default_verbose (void)
{
  return bspline_default_verbose;
}

