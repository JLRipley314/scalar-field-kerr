/*
 * Routines for computing pseudospectral Chebyshev derivatives over
 * interval [lower_, upper_]
 *
 */ 
#include <stddef.h>
#include <math.h>
#include <assert.h>
#include <fftw3.h>
/*
 * global variables
 */
static size_t n_;

static double lower_;
static double upper_;
static double jacobian_; 
static double *pts_; /* Chebyshev collocation
                        points over interval [lower_,upper_]. */

static double *low_pass_; /* For low pass filter in Chebyshev space. */

static double *internal_ch_; /* Stand-in for chebyshev coefficients 
                                and position space for some routines. */

static fftw_plan plan_dct_; /* For fftw Fourier transform. */
/*==========================================================================*/
/* setup fftw and vectors */
/*==========================================================================*/
void init(const size_t n, const double lower, const double upper) 
{
   n_ = n;
   lower_ = lower;
   upper_ = upper;
   jacobian_ = (upper_-lower_) / 2.0;

   pts_         = (double *)fftw_malloc(n_*sizeof(double));
   low_pass_    = (double *)fftw_malloc(n_*sizeof(double));
   internal_ch_ = (double *)fftw_malloc(n_*sizeof(double));
   /* 
    * Chebyshev points on interval [lower, upper]
    */
   for (size_t i=0; i<n_; i++) {
      pts_[i] = jacobian_*cos(M_PI*i/(n_-1.0)) + ((upper_+lower_)/2.0);
   }
   for (size_t i=0; i<n_; i++) {
      low_pass_[i] = exp(-40.0*pow((double)i/n_,10));
   }   
   /* 
    * initialize fftw discrete cosine transform plan 
    */
   double *in_  = (double *)fftw_malloc(n_*sizeof(double));
   double *out_ = (double *)fftw_malloc(n_*sizeof(double));
   assert(in_ !=NULL);
   assert(out_!=NULL);

   for (size_t i=0; i<n_; i++) {
      in_[i]  = 0;
      out_[i] = 0;
   }
   plan_dct_ = fftw_plan_r2r_1d(n_,in_,out_,FFTW_REDFT00,FFTW_PATIENT);
   assert(plan_dct_!=NULL);

   fftw_free(in_);
   fftw_free(out_);
}
/*==========================================================================*/
void cleanup()
{
   assert(pts_        !=NULL);
   assert(low_pass_   !=NULL);
   assert(internal_ch_!=NULL);

   fftw_free(pts_);
   fftw_free(low_pass_);
   fftw_free(internal_ch_);

   assert(plan_dct_!=NULL);
   fftw_destroy_plan(plan_dct_);
   fftw_cleanup();
}
/*==========================================================================*/
/* position space to Chebyshev space */
/*==========================================================================*/
void to_ch(double *po, double *ch) 
{
   assert(po!=NULL);
   assert(ch!=NULL);

   fftw_execute_r2r(plan_dct_, po, ch);
   /*
    * normalize Chebyshev coefficients
    */
   ch[0]    /= (2.0*(n_-1));
   ch[n_-1] /= (2.0*(n_-1));
   for (size_t i=1; i<n_-1; i++) {
      ch[i] /= (n_-1.0);
   }
}
/*==========================================================================*/
/* Chebyshev space to position space */
/*==========================================================================*/
void to_po(double *ch, double *po) 
{
   assert(ch!=NULL);
   assert(po!=NULL);

   /*
    * normalize coefficients for Fourier transform 
    */
   for (size_t i=1; i<n_-1; i++) {
      ch[i] /= 2.0;
   }
   fftw_execute_r2r(plan_dct_, ch, po);
}
/*==========================================================================*/
/* Compute derivative over interval */
/*==========================================================================*/
void der(double *v, double *dv)
{
   assert(v!=NULL);
   assert(dv!=NULL);

   to_ch(v,internal_ch_);
   /*
    * to start Chebyshev derivative recurrence relation 
    */
   internal_ch_[n_-1] = 0;
   internal_ch_[n_-2] = 0;
   /* 
    * use dv as a temporary array
    */
   for (size_t i=0; i<n_; i++) {
      dv[i] = internal_ch_[i];
   }
   /* 
    * apply Chebyshev derivative recurrence relation 
    */
   for (size_t i=n_-2; i>=1; i--) { 
      internal_ch_[i-1] = 2.0*i*dv[i] + internal_ch_[i+1];
   } 
   internal_ch_[0] /= 2.0;
   /* 
    * Normalize derivative to inverval 
    */
   to_po(internal_ch_,dv);
   for (size_t i=0; i<n_; i++) {
      dv[i] /= jacobian_;
   }
}
/*==========================================================================*/
/* Low pass filter of Chebyshev coefficients */
/*==========================================================================*/
void filter(double *v)
{
   assert(v!=NULL);

   to_ch(v,internal_ch_);
   for (size_t i=0; i<n_; i++) { 
      internal_ch_[i] *= low_pass_[i];
   } 
   to_po(internal_ch_,v);
}
/*==========================================================================*/
size_t n() { return n_; }
double lower() { return lower_; }
double upper() { return upper_; }
double pt(const size_t i) { return pts_[i]; }
/*==========================================================================*/
