/*

 * ESPS was originally developed by Entropic Inc., which was later acquired by
 * Microsoft in 1999.  Eventually, the source code was released under the BSD
 * license by KTH.

 * Original info and license below:

 * This material contains unpublished, proprietary software of
 * Entropic Research Laboratory, Inc. Any reproduction, distribution,
 * or publication of this work must be authorized in writing by Entropic
 * Research Laboratory, Inc., and must bear the notice:
 *
 *    "Copyright (c) 1990-1996 Entropic Research Laboratory, Inc.
 *                   All rights reserved"
 *
 * The copyright notice above does not evidence any actual or intended
 * publication of this source code.    
 *
 * Written by:  David Talkin
 * Checked by:
 * Revised by:  Derek Lin, David Talkin
 *
 * Brief description:
 *     A collection of pretty generic signal-processing routines.
 *
 *
 */

//static char *sccs_id = "@(#)sigproc.c 1.4     9/9/96  ERL";

#include "stdafx.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "speech_tools.h"
#include "f0.h"
#include "f0_structs.h"

#define TRUE 1
#define FALSE 0

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Return a time-weighting window of type type and length n in dout.
 * Dout is assumed to be at least n elements long.  Type is decoded in
 * the switch statement below.
 */
int get_window(
    get_f0_session *session,
    register float *dout,
    register int n,
    register int type)
{
    float preemp = 0.0;

    if(n > session->sp.n0) {
        register float *p;
        register int i;

        if (session->sp.din) {
            free(session->sp.din);
            session->sp.din = 0;
        }
        session->sp.din = (float*) malloc(sizeof(float) * n);
        if(! session->sp.din) {
            Fprintf(stderr,"Allocation problems in get_window()\n");
            return 0;
        }
        for(i=0, p=session->sp.din; i++ < n; )
            *p++ = 1;
        session->sp.n0 = n;
    }
    return window(session, session->sp.din, dout, n, preemp, type);
}
 
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Apply a rectangular window (i.e. none).  Optionally, preemphasize. */
void rwindow(
        register float *din,
        register float *dout,
        register int n,
        register float preemp
        )
{
    register float *p;

    /* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
    if(preemp != 0.0) {
        for( p=din+1; n-- > 0; )
            *dout++ = (float)(*p++) - (preemp * *din++);
    } else {
        for( ; n-- > 0; )
            *dout++ =  *din++;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generate a cos^4 window, if one does not already exist. */
void cwindow(
    get_f0_session *session,
    register float *din,
    register float *dout,
    register int n,
    register float preemp)
{
    register int i;
    register float *p;
    register float *q, co;

    if (session->sp.c_wsize != n) {             /* Need to create a new cos**4 window? */
        register double arg, half=0.5;

        if (session->sp.c_wind)
            session->sp.c_wind = (float*) realloc(session->sp.c_wind,
                                                  n*sizeof(float));
        else
            session->sp.c_wind = (float*) malloc(n*sizeof(float));
        session->sp.c_wsize = n;
        arg = 3.1415927 * 2.0 / session->sp.c_wsize;
        for(i=0, q=session->sp.c_wind; i < n; ) {
            co = half*(1.0 - cos((half + (double)i++) * arg));
            *q++ = co * co * co * co;
        }
    }
    /* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
    if(preemp != 0.0) {
        for(i=n, p=din+1, q=session->sp.c_wind; i--; )
            *dout++ = *q++ * ((float)(*p++) - (preemp * *din++));
    } else {
        for(i=n, q=session->sp.c_wind; i--; )
            *dout++ = *q++ * *din++;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generate a Hamming window, if one does not already exist. */
void hwindow(
        get_f0_session *session,
        register float *din,
        register float *dout,
        register int n,
        register float preemp
        )
{
    register int i;
    register float *p;
    register float *q;

    if(session->sp.h_wsize != n) {              /* Need to create a new Hamming window? */
        register double arg, half=0.5;

        if (session->sp.h_wind)
            session->sp.h_wind = (float*)realloc(session->sp.h_wind,n*sizeof(float));
        else
            session->sp.h_wind = (float*)malloc(n*sizeof(float));
        session->sp.h_wsize = n;
        for(i=0, arg=3.1415927*2.0/(session->sp.h_wsize), q=session->sp.h_wind; i < n; )
            *q++ = (.54 - .46 * cos((half + (double)i++) * arg));
    }
    /* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
    if(preemp != 0.0) {
        for(i=n, p=din+1, q=session->sp.h_wind; i--; )
            *dout++ = *q++ * ((float)(*p++) - (preemp * *din++));
    } else {
        for(i=n, q=session->sp.h_wind; i--; )
            *dout++ = *q++ * *din++;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generate a Hanning window, if one does not already exist. */
void hnwindow(
        get_f0_session *session,
        register float *din,
        register float *dout,
        register int n,
        register float preemp
        )
{
    register int i;
    register float *p;
    register float *q;

    if (session->sp.hn_wsize != n) {            /* Need to create a new Hanning window? */
        register double arg, half=0.5;

        if (session->sp.hn_wind)
            session->sp.hn_wind = (float*)realloc(session->sp.hn_wind,
                                                  n*sizeof(float));
        else
            session->sp.hn_wind = (float*)malloc(n*sizeof(float));
        session->sp.hn_wsize = n;
        for(i=0, arg=3.1415927*2.0/(session->sp.hn_wsize), q=session->sp.hn_wind; i < n; )
            *q++ = (half - half * cos((half + (double)i++) * arg));
    }
    /* If preemphasis is to be performed,  this assumes that there are n+1 valid
   samples in the input buffer (din). */
    if(preemp != 0.0) {
        for(i=n, p=din+1, q=session->sp.hn_wind; i--; )
            *dout++ = *q++ * ((float)(*p++) - (preemp * *din++));
    } else {
        for(i=n, q=session->sp.hn_wind; i--; )
            *dout++ = *q++ * *din++;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Apply a window of type type to the short PCM sequence of length n
 * in din.  Return the floating-point result sequence in dout.  If preemp
 * is non-zero, apply preemphasis to tha data as it is windowed.
 */
int window(
    get_f0_session *session,
    register float *din,
    register float *dout,
    register int n,
    register float preemp,
    int type)
{
    switch(type) {
    case 0:                     /* rectangular */
        rwindow(din, dout, n, preemp);
        break;
    case 1:                     /* Hamming */
        hwindow(session, din, dout, n, preemp);
        break;
    case 2:                     /* cos^4 */
        cwindow(session, din, dout, n, preemp);
        break;
    case 3:                     /* Hanning */
        hnwindow(session, din, dout, n, preemp);
        break;
    default:
        Fprintf(stderr,"Unknown window type (%d) requested in window()\n",type);
        return(FALSE);
    }
    return(TRUE);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Compute the pp+1 autocorrelation lags of the windowsize samples in s.
 * Return the normalized autocorrelation coefficients in r.
 * The rms is returned in e.
 */
void autoc(
        register int windowsize,
        register float *s,
        register int p,
        register float *r,
        register float *e
        )
{
    register int i, j;
    register float *q, *t, sum, sum0;

    for( i=windowsize, q=s, sum0=0.0; i--;) {
        sum = *q++;
        sum0 += sum*sum;
    }
    *r = 1.;                    /* r[0] will always =1. */
    if(sum0 == 0.0) {           /* No energy: fake low-energy white noise. */
        *e = 1.;                        /* Arbitrarily assign 1 to rms. */
        /* Now fake autocorrelation of white noise. */
        for ( i=1; i<=p; i++){
            r[i] = 0.;
        }
        return;
    }
    *e = sqrt((double)(sum0/windowsize));
    sum0 = 1.0/sum0;
    for( i=1; i <= p; i++){
        for( sum=0.0, j=windowsize-i, q=s, t=s+i; j--; )
            sum += (*q++) * (*t++);
        *(++r) = sum*sum0;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Using Durbin's recursion, convert the autocorrelation sequence in r
 * to reflection coefficients in k and predictor coefficients in a.
 * The prediction error energy (gain) is left in *ex.
 * Note: durbin returns the coefficients in normal sign format.
 *      (i.e. a[0] is assumed to be = +1.)
 */
void durbin (
        register float *r,
        register float *k,
        register float *a,
        register int p,    /* analysis order */
        register float *ex
        )
{
    float  bb[BIGSORD];
    register int i, j;
    register float e, s, *b = bb;

    e = *r;
    *k = -r[1]/e;
    *a = *k;
    e *= (1. - (*k) * (*k));
    for ( i=1; i < p; i++){
        s = 0;
        for ( j=0; j<i; j++){
            s -= a[j] * r[i-j];
        }
        k[i] = ( s - r[i+1] )/e;
        a[i] = k[i];
        for ( j=0; j<=i; j++){
            b[j] = a[j];
        }
        for ( j=0; j<i; j++){
            a[j] += k[i] * b[i-j-1];
        }
        e *= ( 1. - (k[i] * k[i]) );
    }
    *ex = e;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*  Compute the autocorrelations of the p LP coefficients in a.
 *  (a[0] is assumed to be = 1 and not explicitely accessed.)
 *  The magnitude of a is returned in c.
 *  2* the other autocorrelation coefficients are returned in b.
 */
void a_to_aca (
        float *a,
        float *b,
        float *c,
        register int p
        )
{
    register float  s, *ap, *a0;
    register int  i, j;

    for ( s=1., ap=a, i = p; i--; ap++ )
        s += *ap * *ap;

    *c = s;
    for ( i = 1; i <= p; i++){
        s = a[i-1];
        for (a0 = a, ap = a+i, j = p-i; j--; )
            s += (*a0++ * *ap++);
        *b++ = 2. * s;
    }
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Compute the Itakura LPC distance between the model represented
 * by the signal autocorrelation (r) and its residual (gain) and
 * the model represented by an LPC autocorrelation (c, b).
 * Both models are of order p.
 * r is assumed normalized and r[0]=1 is not explicitely accessed.
 * Values returned by the function are >= 1.
 */
float itakura (
        register int p,
        register float *b,
        register float *c,
        register float *r,
        register float *gain
        )
{
    register float s;

    for( s= *c; p--; )
        s += *r++ * *b++;

    return (s/ *gain);
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Compute the time-weighted RMS of a size segment of data.  The data
 * is weighted by a window of type w_type before RMS computation.  w_type
 * is decoded above in window().
 */
float wind_energy(
    get_f0_session *session,
    register float *data,  /* input PCM data */
    register int size,     /* size of window */
    register int w_type    /* window type */
    )
{
    register float *dp, sum, f;
    register int i;

    if (session->sp.we_nwind < size) {
        if (session->sp.we_dwind)
            session->sp.we_dwind = (float*)realloc(session->sp.we_dwind,
                                                   size*sizeof(float));
        else
            session->sp.we_dwind = (float*)malloc(size*sizeof(float));
        if (!session->sp.we_dwind) {
            Fprintf(stderr,"Can't allocate scratch memory in wind_energy()\n");
            return(0.0);
        }
    }
    if(session->sp.we_nwind != size) {
        get_window(session, session->sp.we_dwind, size, w_type);
        session->sp.we_nwind = size;
    }

    for(i=size, dp = session->sp.we_dwind, sum = 0.0; i-- > 0; ) {
        f = *dp++ * (float)(*data++);
        sum += f*f;
    }
    return((float)sqrt((double)(sum/size)));
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Generic autocorrelation LPC analysis of the short-integer data
 * sequence in data.
 */
int lpc(
    get_f0_session *session,
    int lpc_ord,      /* Analysis order */
    float lpc_stabl,  /* Stability factor to prevent numerical problems. */
    int wsize,        /* window size in points */
    float *data,      /* input data sequence; assumed to be wsize+1 long */
    float *lpca,      /* if non-NULL, return vvector for predictors */
    float *ar,        /* if non-NULL, return vector for normalized autoc. */
    float *lpck,      /* if non-NULL, return vector for PARCOR's */
    float *normerr,   /* return scaler for normalized error */
    float *rms,       /* return scaler for energy in preemphasized window */
    float preemp,
    int type          /* window type (decoded in window() above) */
    )
{
    float rho[BIGSORD+1], k[BIGSORD], a[BIGSORD+1],*r,*kp,*ap,en,er,wfact;

    if((wsize <= 0) || (!data) || (lpc_ord > BIGSORD)) return(FALSE);

    if (session->sp.lpc_nwind != wsize) {
        if(session->sp.lpc_dwind)
            session->sp.lpc_dwind = (float*)realloc(session->sp.lpc_dwind,
                                                    wsize*sizeof(float));
        else
            session->sp.lpc_dwind = (float*)malloc(wsize*sizeof(float));
        if (!session->sp.lpc_dwind) {
            Fprintf(stderr,"Can't allocate scratch memory in lpc()\n");
            return(FALSE);
        }
        session->sp.lpc_nwind = wsize;
    }

    window(session, data, session->sp.lpc_dwind, wsize, preemp, type);
    if(!(r = ar)) r = rho;      /* Permit optional return of the various */
    if(!(kp = lpck)) kp = k;    /* coefficients and intermediate results. */
    if(!(ap = lpca)) ap = a;
    autoc( wsize, session->sp.lpc_dwind, lpc_ord, r, &en );
    if(lpc_stabl > 1.0) {       /* add a little to the diagonal for stability */
        int i;
        float ffact;
        ffact =1.0/(1.0 + exp((-lpc_stabl/20.0) * log(10.0)));
        for(i=1; i <= lpc_ord; i++) rho[i] = ffact * r[i];
        *rho = *r;
        r = rho;
        if(ar)
            for(i=0;i<=lpc_ord; i++) ar[i] = r[i];
    }
    durbin ( r, kp, &ap[1], lpc_ord, &er);
    switch(type) {              /* rms correction for window */
    case 0:
        wfact = 1.0;            /* rectangular */
        break;
    case 1:
        wfact = .630397;                /* Hamming */
        break;
    case 2:
        wfact = .443149;                /* (.5 - .5*cos)^4 */
        break;
    case 3:
        wfact = .612372;                /* Hanning */
        break;
    default:
        wfact = 1.0;
        break;
    }
    *ap = 1.0;
    if(rms) *rms = en/wfact;
    if(normerr) *normerr = er;
    return(TRUE);
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Return a sequence based on the normalized crosscorrelation of the signal
   in data.
 *
  data is the input speech array
  size is the number of samples in each correlation
  start is the first lag to compute (governed by the highest expected F0)
  nlags is the number of cross correlations to compute (set by lowest F0)
  engref is the energy computed at lag=0 (i.e. energy in ref. window)
  maxloc is the lag at which the maximum in the correlation was found
  maxval is the value of the maximum in the CCF over the requested lag interval
  correl is the array of nlags cross-correlation coefficients (-1.0 to 1.0)
 *
 */
int crossf(
        get_f0_session *session,
        float *data,
        int size,
        int start,
        int nlags,
        float *engref,
        int *maxloc,
        float *maxval,
        float *correl
        )
{
    register float *dp, *ds, sum, st;
    register int j;
    register  float *dq, t, *p, engr, *dds, amax;
    register  double engc;
    int i, iloc, total;
    //int sizei, sizeo, maxsize;

    /* Compute mean in reference window and subtract this from the
     entire sequence.  This doesn't do too much damage to the data
     sequenced for the purposes of F0 estimation and removes the need for
     more principled (and costly) low-cut filtering. */
    if((total = size+start+nlags) > session->sp.f_dbsize) {
        if(session->sp.f_dbdata)
            free(session->sp.f_dbdata);
        session->sp.f_dbdata = 0;
        session->sp.f_dbsize = 0;
        session->sp.f_dbdata = (float*)malloc(sizeof(float)*total);
        if (!session->sp.f_dbdata) {
            Fprintf(stderr,"Allocation failure in crossf()\n");
            return 1;
        }
        session->sp.f_dbsize = total;
    }
    for(engr=0.0, j=size, p=data; j--; )
        engr += *p++;
    engr /= size;
    for(j=size+nlags+start, dq = session->sp.f_dbdata, p=data; j--; )
        *dq++ = *p++ - engr;

    //maxsize = start + nlags;
    //sizei = size + start + nlags + 1;
    //sizeo = nlags + 1;

    /* Compute energy in reference window. */
    for(j=size, dp=session->sp.f_dbdata, sum=0.0; j--; ) {
        st = *dp++;
        sum += st * st;
    }

    *engref = engr = sum;
    if(engr > 0.0) {    /* If there is any signal energy to work with... */
        /* Compute energy at the first requested lag. */
        for(j=size, dp=session->sp.f_dbdata+start, sum=0.0; j--; ) {
            st = *dp++;
            sum += st * st;
        }
        engc = sum;

        /* COMPUTE CORRELATIONS AT ALL OTHER REQUESTED LAGS. */
        for(i=0, dq=correl, amax=0.0, iloc = -1; i < nlags; i++) {
            dp = session->sp.f_dbdata;
            dds = ds = session->sp.f_dbdata + i + start;
            for(j=size, sum=0.0; j--; )
                sum += *dp++ * *ds++;
            *dq++ = t = (sum/sqrt((double)(engc*engr))); /* output norm. CC */
            engc -= (double)(*dds * *dds); /* adjust norm. energy for next lag */
            if((engc += (double)(*ds * *ds)) < 1.0)
                engc = 1.0;             /* (hack: in case of roundoff error) */
            if(t > amax) {              /* Find abs. max. as we go. */
                amax = t;
                iloc = i+start;
            }
        }
        *maxloc = iloc;
        *maxval = amax;
    } else {    /* No energy in signal; fake reasonable return vals. */
        *maxloc = 0;
        *maxval = 0.0;
        for(p=correl,i=nlags; i-- > 0; )
            *p++ = 0.0;
    }

    return 0;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/* Return a sequence based on the normalized crosscorrelation of the
   signal in data.  This is similar to crossf(), but is designed to
   compute only small patches of the correlation sequence.  The length of
   each patch is determined by nlags; the number of patches by nlocs, and
   the locations of the patches is specified by the array locs.  Regions
   of the CCF that are not computed are set to 0.
 *
  data is the input speech array
  size is the number of samples in each correlation
  start0 is the first (virtual) lag to compute (governed by highest F0)
  nlags0 is the number of lags (virtual+actual) in the correlation sequence
  nlags is the number of cross correlations to compute at each location
  engref is the energy computed at lag=0 (i.e. energy in ref. window)
  maxloc is the lag at which the maximum in the correlation was found
  maxval is the value of the maximum in the CCF over the requested lag interval
  correl is the array of nlags cross-correlation coefficients (-1.0 to 1.0)
  locs is an array of indices pointing to the center of a patches where the
       cross correlation is to be computed.
  nlocs is the number of correlation patches to compute.
 *
 */
int crossfi(
        get_f0_session *session,
        float *data,
        int size,
        int start0,
        int nlags0,
        int nlags,
        float *engref,
        int *maxloc,
        float *maxval,
        float *correl,
        int *locs,
        int nlocs
        )
{
    register float *dp, *ds, sum, st;
    register int j;
    register  float *dq, t, *p, engr, *dds, amax;
    register  double engc;
    int i, iloc, start, total;

    /* Compute mean in reference window and subtract this from the
     entire sequence. */
    if((total = size+start0+nlags0) > session->sp.fi_dbsize) {
        if(session->sp.fi_dbdata)
            free(session->sp.fi_dbdata);
        session->sp.fi_dbdata = 0;
        session->sp.fi_dbsize = 0;
        if (!(session->sp.fi_dbdata = (float*)malloc(sizeof(float)*total))) {
            Fprintf(stderr,"Allocation failure in crossfi()\n");
            return 1;
        }
        session->sp.fi_dbsize = total;
    }
    for(engr=0.0, j=size, p=data; j--; ) engr += *p++;
    engr /= size;
    for(j=size+nlags0+start0, dq = session->sp.fi_dbdata, p=data; j--; ) {
        *dq++ = *p++ - engr;
    }

    /* Zero the correlation output array to avoid confusing the peak
     picker (since all lags will not be computed). */
    for(p=correl,i=nlags0; i-- > 0; )
        *p++ = 0.0;

    /* compute energy in reference window */
    for(j=size, dp=session->sp.fi_dbdata, sum=0.0; j--; ) {
        st = *dp++;
        sum += st * st;
    }

    *engref = engr = sum;
    amax=0.0;
    iloc = -1;
    if(engr > 0.0) {
        for( ; nlocs > 0; nlocs--, locs++ ) {
            start = *locs - (nlags>>1);
            if(start < start0)
                start = start0;
            dq = correl + start - start0;
            /* compute energy at first requested lag */
            for(j=size, dp=session->sp.fi_dbdata+start, sum=0.0; j--; ) {
                st = *dp++;
                sum += st * st;
            }
            engc = sum;

            /* COMPUTE CORRELATIONS AT ALL REQUESTED LAGS */
            for(i=0; i < nlags; i++) {
                dp=session->sp.fi_dbdata;
                dds = ds = session->sp.fi_dbdata + i + start;
                for(j=size, sum=0.0; j--; )
                    sum += *dp++ * *ds++;
                if(engc < 1.0)
                    engc = 1.0;         /* in case of roundoff error */
                *dq++ = t = (sum/sqrt((double)(10000.0 + (engc*engr))));
                engc -= (double)(*dds * *dds);
                engc += (double)(*ds * *ds);
                if(t > amax) {
                    amax = t;
                    iloc = i+start;
                }
            }
        }
        *maxloc = iloc;
        *maxval = amax;
    } else {
        *maxloc = 0;
        *maxval = 0.0;
    }
    return 0;
}

