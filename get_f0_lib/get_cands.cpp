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
 * Revised by:
 *
 * Brief description:
 *
 */

//static char *sccs_id = "@(#)get_cands.c       1.5     9/9/96  ERL";

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>

#include "speech_tools.h"
#include "f0_structs.h"
#include "f0.h"

#define TRUE 1
#define FALSE 0

static void get_cand(Cross *cross, float *peak, int *loc, int nlags, int *ncand, float cand_thresh);
static void peak(float *y, float *xp, float *yp);
static int do_ffir(
    get_f0_session *session,
    register float *buf,
    register int in_samps,
    register int *out_samps,
    int idx,
    register int invert,
    register int skip,
    register int init);
static int lc_lin_fir(register float fc, int *nf, float *coef);
static int downsamp(
    get_f0_session *session,
    float *in,
    int samples,
    int *outsamps,
    int state_idx,
    int decimate,
    int init);

/* ----------------------------------------------------------------------- */
int get_fast_cands(
        get_f0_session *session,
        float *fdata,
        float *fdsdata,
        int ind,
        int step,
        int size,
        int dec,
        int start,
        int nlags,
        float *engref,
        int *maxloc,
        float *maxval,
        Cross *cp,
        float *peaks,
        int *locs,
        int *ncand
        )
{
    int decind, decstart, decnlags, decsize, i, j, *lp;
    float *corp, xp, yp, lag_wt;
    register float *pe;
    F0_params *par = session->par;

    lag_wt = par->lag_weight/nlags;
    decnlags = 1 + (nlags/dec);
    if((decstart = start/dec) < 1) decstart = 1;
    decind = (ind * step)/dec;
    decsize = 1 + (size/dec);
    corp = cp->correl;
   

    if (crossf(session, fdsdata + decind, decsize, decstart, decnlags,
               engref, maxloc, maxval, corp) != 0) {
        return 1; // failure
    }

    cp->maxloc = *maxloc;       /* location of maximum in correlation */
    cp->maxval = *maxval;       /* max. correlation value (found at maxloc) */
    cp->rms = sqrt(*engref/size); /* rms in reference window */
    cp->firstlag = decstart;

    get_cand(cp,peaks,locs,decnlags,ncand,par->cand_thresh); /* return high peaks in xcorr */

    /* Interpolate to estimate peak locations and values at high sample rate. */
    for(i = *ncand, lp = locs, pe = peaks; i--; pe++, lp++) {
        j = *lp - decstart - 1;
        peak(&corp[j],&xp,&yp);
        *lp = (*lp * dec) + (int)(0.5+(xp*dec)); /* refined lag */
        *pe = yp*(1.0 - (lag_wt* *lp)); /* refined amplitude */
    }

    if(*ncand >= par->n_cands) {        /* need to prune candidates? */
        register int *loc, *locm, lt;
        register float smaxval, *pem;
        register int outer, inner, lim;
        for(outer=0, lim = par->n_cands-1; outer < lim; outer++)
            for(inner = *ncand - 1 - outer,
                pe = peaks + (*ncand) -1, pem = pe-1,
                loc = locs + (*ncand) - 1, locm = loc-1;
                inner--;
                pe--,pem--,loc--,locm--)
                if((smaxval = *pe) > *pem) {
                    *pe = *pem;
                    *pem = smaxval;
                    lt = *loc;
                    *loc = *locm;
                    *locm = lt;
                }
        *ncand = par->n_cands-1;  /* leave room for the unvoiced hypothesis */
    }

    if (crossfi(session, fdata + (ind * step), size, start, nlags, 7,
                engref, maxloc, maxval, corp, locs, *ncand) != 0) {
        return 1; // failure
    }

    cp->maxloc = *maxloc;       /* location of maximum in correlation */
    cp->maxval = *maxval;       /* max. correlation value (found at maxloc) */
    cp->rms = sqrt(*engref/size); /* rms in reference window */
    cp->firstlag = start;
    get_cand(cp,peaks,locs,nlags,ncand,par->cand_thresh); /* return high peaks in xcorr */
    if(*ncand >= par->n_cands) {        /* need to prune candidates again? */
        register int *loc, *locm, lt;
        register float smaxval, *pe, *pem;
        register int outer, inner, lim;
        for(outer=0, lim = par->n_cands-1; outer < lim; outer++)
            for(inner = *ncand - 1 - outer,
                pe = peaks + (*ncand) -1, pem = pe-1,
                loc = locs + (*ncand) - 1, locm = loc-1;
                inner--;
                pe--,pem--,loc--,locm--)
                if((smaxval = *pe) > *pem) {
                    *pe = *pem;
                    *pem = smaxval;
                    lt = *loc;
                    *loc = *locm;
                    *locm = lt;
                }
        *ncand = par->n_cands - 1;  /* leave room for the unvoiced hypothesis */
    }
    return 0; // success
}

/* ----------------------------------------------------------------------- */
float *downsample(
    get_f0_session *session,
    float *input,
    int samsin,
    int state_idx,
    double freq,
    int *samsout,
    int decimate,
    int last_time
    )
{
    float       beta = 0.0;
    int init;

    if(input && (samsin > 0) && (decimate > 0) && *samsout) {
        if(decimate == 1) {
            return(input);
        }

        if (session->dp.first_time){
            int nbuff = (samsin/decimate) + (2*session->gcs.ncoeff);

            session->gcs.ncoeff = ((int)(freq * .005)) | 1;
            beta = .5/decimate;
            session->gcs.foutput = (float*)malloc(sizeof(float) * nbuff);
            if (!session->gcs.foutput) {
                fprintf(stderr, "Can't allocate foutput in downsample");
                return 0;
            }

            for( ; nbuff > 0 ;)
                session->gcs.foutput[--nbuff] = 0.0;

            if (!lc_lin_fir(beta, &session->gcs.ncoeff, session->gcs.b)) {
                fprintf(stderr,"\nProblems computing interpolation filter\n");
                free(session->gcs.foutput);
                session->gcs.foutput = 0;
                return 0;
            }
            session->gcs.ncoefft = (session->gcs.ncoeff/2) + 1;
        }                   /*  endif new coefficients need to be computed */

        if (session->dp.first_time)
            init = 1;
        else if (last_time)
            init = 2;
        else
            init = 0;

        if (downsamp(
                    session,
                    input,
                    samsin,
                    samsout,
                    state_idx,
                    decimate,
                    init
                    ) == TRUE)
            return session->gcs.foutput;
        else
            Fprintf(stderr,"Problems in downsamp() in downsample()\n");
    } else
        Fprintf(stderr,"Bad parameters passed to downsample()\n");

    return(NULL);
}

/* ----------------------------------------------------------------------- */
/* Get likely candidates for F0 peaks. */
static void get_cand(
                Cross *cross,
                float *peak,
                int *loc,
                int nlags,
                int *ncand,
                float cand_thresh)
{
  register int i, lastl, *t;
  register float o, p, q, *r, *s, clip;
  int start, ncan; //, maxl;

  clip = cand_thresh * cross->maxval;
  //maxl = cross->maxloc;
  lastl = nlags - 2;
  start = cross->firstlag;

  r = cross->correl;
  o= *r++;                      /* first point */
  q = *r++;                     /* middle point */
  p = *r++;
  s = peak;
  t = loc;
  ncan=0;
  for(i=1; i < lastl; i++, o=q, q=p, p= *r++){
    if((q > clip) &&            /* is this a high enough value? */
      (q >= p) && (q >= o)){ /* NOTE: this finds SHOLDERS and PLATEAUS
                                      as well as peaks (is this a good idea?) */
        *s++ = q;               /* record the peak value */
        *t++ = i + start;       /* and its location */
        ncan++;                 /* count number of peaks found */
      }
  }
/*
  o = q;
  q = p;
  if( (q > clip) && (q >=0)){
    *s++ = q;
    *t++ = i+start;
    ncan++;
  }
*/
  *ncand = ncan;
}

/* ----------------------------------------------------------------------- */
/* buffer-to-buffer downsample operation */
/* This is STRICTLY a decimator! (no upsample) */
static int downsamp(
    get_f0_session *session,
    float *in,
    int samples,
    int *outsamps,
    int state_idx,
    int decimate,
    int init)
{
    if(in && session->gcs.foutput) {
        if (do_ffir(session, in, samples,
                       outsamps, state_idx, 0,
                       decimate, init) == 0)
            return TRUE;
        else
            return FALSE;
    } else
        fprintf(stderr, "Bad signal(s) passed to downsamp()\n");
    return FALSE;
}

/*      ----------------------------------------------------------      */
static int do_ffir(
    get_f0_session *session,
    register float *buf,
    register int in_samps,
    register int *out_samps,
    int idx,
    register int invert,
    register int skip,
    register int init)
/* fc contains 1/2 the coefficients of a symmetric FIR filter with unity
    passband gain.  This filter is convolved with the signal in buf.
    The output is placed in buf2.  If(invert), the filter magnitude
    response will be inverted.  If(init&1), beginning of signal is in buf;
    if(init&2), end of signal is in buf.  out_samps is set to the number of
    output points placed in bufo. */
{
    register float *dp1, *dp2, *dp3, sum, integral;
    register int i, j, k, l;
    register float *sp;
    register float *buf1;
    register float *bufo = session->gcs.foutput;

    buf1 = buf;
    if(session->gcs.ncoefft > session->gcs.fsize) {/*allocate memory for full coeff. array and filter memory */
        if (session->gcs.co) {
            free(session->gcs.co);
            session->gcs.co = 0;
        }
        if (session->gcs.mem) {
            free(session->gcs.mem);
            session->gcs.mem = 0;
        }
        session->gcs.fsize = 0;
        i = (session->gcs.ncoefft+1) * 2;

        session->gcs.co = (float *) malloc(sizeof(float)*i);
        if (session->gcs.co == 0) {
            fprintf(stderr,"allocation problems in do_fir()\n");
            return 1; // error
        }
        session->gcs.mem = (float *) malloc(sizeof(float)*i);
        if (session->gcs.mem == 0) {
            fprintf(stderr,"allocation problems in do_fir()\n");
            free(session->gcs.co);
            session->gcs.co = 0;
            return 1; // error
        }

        session->gcs.fsize = session->gcs.ncoefft;
    }

    /* fill 2nd half with data */
    i = session->gcs.ncoefft;
    dp1 = session->gcs.mem + session->gcs.ncoefft - 1;
    for(; i-- > 0; )  *dp1++ = *buf++;

    if(init & 1) {      /* Is the beginning of the signal in buf? */
        /* Copy the half-filter and its mirror image into the coefficient array. */
        i = session->gcs.ncoefft - 1;
        dp1 = session->gcs.co + (session->gcs.ncoefft - 1) * 2;
        dp2 = session->gcs.co;
        dp3 = session->gcs.b + session->gcs.ncoefft - 1;
        for (integral = 0.0; i-- > 0; ) {
            if (!invert)
                *dp1-- = *dp2++ = *dp3--;
            else {
                integral += (sum = *dp3--);
                *dp1-- = *dp2++ = -sum;
            }
        }

        if (!invert)
            *dp1 = *dp3;        /* point of symmetry */
        else {
            integral *= 2;
            integral += *dp3;
            *dp1 = integral - *dp3;
        }

        for (i=session->gcs.ncoefft-1, dp1=session->gcs.mem; i-- > 0; )
            *dp1++ = 0;
    }
    else {
        i = session->gcs.ncoefft - 1;
        dp1 = session->gcs.mem;
        sp = session->gcs.state;
        for(; i-- > 0; )
            *dp1++ = *sp++;
    }

    i = in_samps;
    session->gcs.resid = 0;

    k = (session->gcs.ncoefft << 1) -1; /* inner-product loop limit */

    if (skip > 1) {  /* skip points (e.g. for downsampling) */
        /* the buffer end is padded with (ncoef-1) data points */
        for( l=0 ; l < *out_samps; l++ ) {
            dp1 = session->gcs.mem;
            dp2 = session->gcs.co;
            dp3 = session->gcs.mem + skip;
            for (j=k-skip, sum=0.0; j-- >0; *dp1++ = *dp3++)
                sum += *dp2++ * *dp1;
            for (j=skip; j-- >0; *dp1++ = *buf++) /* new data to memory */
                sum += *dp2++ * *dp1;
            *bufo++ = (sum<0.0) ? sum -0.5 : sum +0.5;
        }
        if(init & 2){
            session->gcs.resid = in_samps - *out_samps * skip;
            for (l=session->gcs.resid/skip; l-- >0; ){
                dp1 = session->gcs.mem;
                dp2 = session->gcs.co;
                dp3 = session->gcs.mem + skip;
                for (j=k-skip, sum=0.0; j-- >0; *dp1++ = *dp3++)
                    sum += *dp2++ * *dp1;
                for (j=skip; j-- >0; *dp1++ = 0.0)
                    sum += *dp2++ * *dp1;
                *bufo++ = (sum<0.0) ? sum -0.5 : sum +0.5;
                (*out_samps)++;
            }
        }
        else {
            dp3 = buf1 + idx - session->gcs.ncoefft + 1;
            l = session->gcs.ncoefft - 1;
            sp = session->gcs.state;
            for(; l-- >0; ) *sp++ = *dp3++;
        }
    }
    return 0;
}

/*      ----------------------------------------------------------      */
static int lc_lin_fir(register float fc, int *nf, float *coef)
/* create the coefficients for a symmetric FIR lowpass filter using the
   window technique with a Hanning window. */
{
    register int        i, n;
    register double     twopi, fn, c;

    if(((*nf % 2) != 1))
        *nf = *nf + 1;
    n = (*nf + 1)/2;

    /*  Compute part of the ideal impulse response (the sin(x)/x kernel). */
    twopi = M_PI * 2.0;
    coef[0] = 2.0 * fc;
    c = M_PI;
    fn = twopi * fc;
    for(i=1;i < n; i++) coef[i] = sin(i * fn)/(c * i);

    /* Now apply a Hanning window to the (infinite) impulse response. */
    /* (Probably should use a better window, like Kaiser...) */
    fn = twopi/(double)(*nf);
    for(i=0;i<n;i++)
        coef[n-i-1] *= (.5 - (.5 * cos(fn * ((double)i + 0.5))));
   
    return(TRUE);
}


/* ----------------------------------------------------------------------- */
/* Use parabolic interpolation over the three points defining the peak
 * vicinity to estimate the "true" peak. */
static void peak(
    float *y,   /* vector of length 3 defining peak */
    float *xp,
    float *yp   /* x,y values of parabolic peak fitting the input points. */
    )
{
    register float a, c;

    a = (y[2]-y[1])+(.5*(y[0]-y[2]));
    if(fabs(a) > .000001) {
        *xp = c = (y[0]-y[2])/(4.0*a);
        *yp = y[1] - (a*c*c);
    } else {
        *xp = 0.0;
        *yp = y[1];
    }
}

