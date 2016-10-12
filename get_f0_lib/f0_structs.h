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
 * @(#)f0_structs.h     1.7 9/9/96 ERL
 * Brief description:
 *
 */


/* f0_structs.h */

#ifndef _F0_STRUCTS_H_
#define _F0_STRUCTS_H_

#define BIGSORD 100

typedef struct cross_rec { /* for storing the crosscorrelation information */
        float   rms;    /* rms energy in the reference window */
        float   maxval; /* max in the crosscorr. fun. q15 */
        short   maxloc; /* lag # at which max occured   */
        short   firstlag; /* the first non-zero lag computed */
        float   *correl; /* the normalized corsscor. fun. q15 */
} Cross;

typedef struct dp_rec { /* for storing the DP information */
        short   ncands; /* # of candidate pitch intervals in the frame */
        short   *locs; /* locations of the candidates */
        float   *pvals; /* peak values of the candidates */
        float   *mpvals; /* modified peak values of the candidates */
        short   *prept; /* pointers to best previous cands. */
        float   *dpvals; /* cumulative error for each candidate */
} Dprec;

typedef struct windstat_rec {  /* for lpc stat measure in a window */
    float rho[BIGSORD+1];
    float err;
    float rms;
} Windstat;

typedef struct sta_rec {  /* for stationarity measure */
  float *stat;
  float *rms;
  float *rms_ratio;
} Stat;


typedef struct frame_rec{
  Cross *cp;
  Dprec *dp;
  float rms;
  struct frame_rec *next;
  struct frame_rec *prev;
} Frame_getf0;

extern Frame_getf0 *alloc_frame(int nlags, int ncands);

#endif
