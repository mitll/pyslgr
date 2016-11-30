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
 * Brief description:  Estimate speech fundamental frequency.
 *
 */

#include "stdafx.h"
#include <stdlib.h>
#include "speech_tools.h"

//static char *sccs_id = "@(#)dp_f0.c   1.14    10/21/96        ERL";

/* A fundamental frequency estimation algorithm using the normalized
   cross correlation function and dynamic programming.  The algorithm
   implemented here is similar to that presented by B. Secrest and
   G. Doddington, "An integrated pitch tracking algorithm for speech
   systems", Proc. ICASSP-83, pp.1352-1355.  It is fully described
   by D. Talkin, "A robust algorithm for ptich tracking (RAPT)", in
   W. B. Kleijn & K. K. Paliwal (eds.) Speech Coding and Synthesis,
   (New York: Elsevier, 1995). */

/* For each frame, up to par->n_cands cross correlation peaks are
   considered as F0 intervals.  Each is scored according to its within-
   frame properties (relative amplitude, relative location), and
   according to its connectivity with each of the candidates in the
   previous frame.  An unvoiced hypothesis is also generated at each
   frame and is considered in the light of voicing state change cost,
   the quality of the cross correlation peak, and frequency continuity. */

/* At each frame, each candidate has associated with it the following
   items:
        its peak value
        its peak value modified by its within-frame properties
        its location
        the candidate # in the previous frame yielding the min. err.
                (this is the optimum path pointer!)
        its cumulative cost: (local cost + connectivity cost +
                cumulative cost of its best-previous-frame-match). */

/* Dynamic programming is then used to pick the best F0 trajectory and voicing
   state given the local and transition costs for the entire utterance. */

/* To avoid the necessity of computing the full crosscorrelation at
   the input sample rate, the signal is downsampled; a full ccf is
   computed at the lower frequency; interpolation is used to estimate the
   location of the peaks at the higher sample rate; and the fine-grained
   ccf is computed only in the vicinity of these estimated peak
   locations. */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include "f0.h"
#include "f0_structs.h"

static int  debug_level = 0;
static char *ProgName = (char *) "get_f0";

/*
 * READ_SIZE: length of input data frame in sec to read
 * DP_CIRCULAR: determines the initial size of DP circular buffer in sec
 * DP_HIST: stored frame history in second before checking for common path
 *      DP_CIRCULAR > READ_SIZE, DP_CIRCULAR at least 2 times of DP_HIST
 * DP_LIMIT: in case no convergence is found, DP frames of DP_LIMIT secs
 *      are kept before output is forced by simply picking the lowest cost
 *      path
 */

#define READ_SIZE 0.2
#define DP_CIRCULAR 1.5
#define DP_HIST 0.5
#define DP_LIMIT 1.0

extern float *downsample(get_f0_session*,
                         float *input,
                         int samsin,
                         int state_idx,
                         double freq,
                         int *samsout,
                         int decimate,
                         int last_time);
static Stat* get_stationarity(get_f0_session*,
                              float *fdata,
                              double freq,
                              int buff_size,
                              int nframes,
                              int frame_step,
                              int first_time);
extern int get_fast_cands(
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
        int *ncand);

extern float wind_energy(
    get_f0_session *session,
    register float *data,
    register int size,
    register int w_type);
extern float itakura (register int p,
                      register float *b,
                      register float *c,
                      register float *r,
                      register float *gain);
void close_dp_f0(get_f0_session*);

/*
 * stationarity parameters -
 * STAT_WSIZE: window size in sec used in measuring frame energy/stationarity
 * STAT_AINT: analysis interval in sec in measuring frame energy/stationarity
 */
#define STAT_WSIZE 0.030
#define STAT_AINT 0.020

/*--------------------------------------------------------------------*/
int get_Nframes(long buffsize, int pad, int step)
{
    if (buffsize < pad)
        return 0;
    else
        return (buffsize - pad) / step;
}


/*--------------------------------------------------------------------*/
int init_dp_f0(
        get_f0_session *session,
        double freq,
        long *buffsize,
        long *sdstep
        )
{
    int nframes;
    int i;
    int stat_wsize, agap, ind, downpatch;
    F0_params *par = session->par;

/*
 * reassigning some constants
 */

  session->dp.tcost = par->trans_cost;
  session->dp.tfact_a = par->trans_amp;
  session->dp.tfact_s = par->trans_spec;
  session->dp.vbias = par->voice_bias;
  session->dp.fdouble = par->double_cost;
  session->dp.frame_int = par->frame_step;

  session->dp.step = Round(session->dp.frame_int * freq);
  session->dp.size = Round(par->wind_dur * freq);
  session->dp.frame_int = ((float)session->dp.step)/freq;
  session->dp.start = Round(freq / par->max_f0);
  session->dp.stop = Round(freq / par->min_f0);
  session->dp.nlags = session->dp.stop - session->dp.start + 1;
  session->dp.ncomp = session->dp.size + session->dp.stop + 1; /* # of samples required by xcorr
                              comp. per fr. */
  //session->dp.maxpeaks = 2 + (session->dp.nlags/2);     /* maximum number of "peaks" findable in ccf */
  // Per the comment in function get_cand, the peak finder accepts shoulders and plateaus
  // So... maxpeaks needs to be bigger
  session->dp.maxpeaks = 2 + session->dp.nlags;  // not sure what the 2 was so I left it...
  session->dp.ln2 = log(2.0);
  session->dp.size_frame_hist = (int) (DP_HIST / session->dp.frame_int);
  session->dp.size_frame_out = (int) (DP_LIMIT / session->dp.frame_int);

/*
 * SET UP THE D.P. WEIGHTING FACTORS:
 *      The intent is to make the effectiveness of the various fudge factors
 *      independent of frame rate or sampling frequency.
 */

  /* Lag-dependent weighting factor to emphasize early peaks (higher freqs)*/
  session->dp.lagwt = par->lag_weight/session->dp.stop;

  /* Penalty for a frequency skip in F0 per frame */
  session->dp.freqwt = par->freq_weight/session->dp.frame_int;

  i = (int) (READ_SIZE *freq);
  if (session->dp.ncomp >= session->dp.step)
    nframes = ((i-session->dp.ncomp)/session->dp.step ) + 1;
  else
    nframes = i / session->dp.step;

  /* *buffsize is the number of samples needed to make F0 computation
     of nframes DP frames possible.  The last DP frame is patched with
     enough points so that F0 computation on it can be carried.  F0
     computaion on each frame needs enough points to do

     1) xcross or cross correlation measure:
           enough points to do xcross - ncomp

     2) stationarity measure:
           enough to make 30 msec windowing possible - ind

     3) downsampling:
           enough to make filtering possible -- downpatch

     So there are nframes whole DP frames, padded with pad points
     to make the last frame F0 computation ok.

  */

  /* last point in data frame needs points of 1/2 downsampler filter length
     long, 0.005 is the filter length used in downsampler */
  downpatch = (((int) (freq * 0.005))+1) / 2;

  stat_wsize = (int) (STAT_WSIZE * freq);
  agap = (int) (STAT_AINT * freq);
  ind = ( agap - stat_wsize ) / 2;
  i = stat_wsize + ind;
  session->dp.pad = downpatch + ((i>session->dp.ncomp) ? i:session->dp.ncomp);
  *buffsize = nframes * session->dp.step + session->dp.pad;
  *sdstep = nframes * session->dp.step;

  /* Allocate space for the DP storage circularly linked data structure */

  session->dp.size_cir_buffer = (int) (DP_CIRCULAR / session->dp.frame_int);

  /* creating circularly linked data structures */
  session->dp.tailF = alloc_frame(session->dp.nlags, par->n_cands);
  if (!session->dp.tailF) {
    return 1;
  }
  session->dp.headF = session->dp.tailF;

  /* link them up */
  for(i=1; i<session->dp.size_cir_buffer; i++){
    session->dp.headF->next = alloc_frame(session->dp.nlags, par->n_cands);
    if (!session->dp.headF->next) {
      close_dp_f0(session);
      return 1;
    }
    session->dp.headF->next->prev = session->dp.headF;
    session->dp.headF = session->dp.headF->next;
  }
  session->dp.headF->next = session->dp.tailF;
  session->dp.tailF->prev = session->dp.headF;

  session->dp.headF = session->dp.tailF;

  /* Allocate sscratch array to use during backtrack convergence test. */
  if( ! session->dp.pcands ) {
    session->dp.pcands = (int *) malloc( par->n_cands * sizeof(int));
    if (! session->dp.pcands) {
      fprintf(stderr, "can't allocate pathcands\n");
      close_dp_f0(session);
      return 1;
    }
  }

  /* Allocate arrays to return F0 and related signals. */

  /* Note: remember to compare *vecsize with size_frame_out, because
     size_cir_buffer is not constant */
  session->dp.output_buf_size = session->dp.size_cir_buffer;
  session->dp.rms_speech = (float*)malloc(sizeof(float) * session->dp.output_buf_size);
  if (!session->dp.rms_speech) {
    fprintf(stderr, "rms_speech malloc failed\n");
    close_dp_f0(session);
    return 1;
  }
  session->dp.f0p = (float*)malloc(sizeof(float) * session->dp.output_buf_size);
  if (!session->dp.f0p) {
    fprintf(stderr, "f0p malloc failed\n");
    close_dp_f0(session);
    return 1;
  }
  session->dp.vuvp = (float*)malloc(sizeof(float)* session->dp.output_buf_size);
  if (!session->dp.vuvp) {
    fprintf(stderr, "vuvp malloc failed\n");
    close_dp_f0(session);
    return 1;
  }
  session->dp.acpkp = (float*)malloc(sizeof(float) * session->dp.output_buf_size);
  if (!session->dp.acpkp) {
    fprintf(stderr, "acpkp malloc failed\n");
    close_dp_f0(session);
    return 1;
  }

  /* Allocate space for peak location and amplitude scratch arrays. */
  session->dp.peaks = (float*)malloc(sizeof(float) * session->dp.maxpeaks);
  if (!session->dp.peaks) {
    fprintf(stderr, "peaks malloc failed\n");
    close_dp_f0(session);
    return 1;
  }
  session->dp.locs = (int*)malloc(sizeof(int) * session->dp.maxpeaks);
  if (!session->dp.locs) {
    fprintf(stderr, "locs malloc failed\n");
    close_dp_f0(session);
    return 1;
  }

  /* Initialise the retrieval/saving scheme of window statistic measures */
  session->dp.wReuse = agap / session->dp.step;
  if (session->dp.wReuse){
      session->dp.windstat = (Windstat *) malloc( session->dp.wReuse * sizeof(Windstat));
      if (!session->dp.windstat) {
        fprintf(stderr, "windstat malloc failed\n");
        close_dp_f0(session);
        return 1;
      }
      for(i=0; i<session->dp.wReuse; i++){
          session->dp.windstat[i].err = 0;
          session->dp.windstat[i].rms = 0;
      }
  }

  if(debug_level){
    Fprintf(stderr, "%s: done with initialization:\n", ProgName);
    Fprintf(stderr,
            " size_cir_buffer:%d  xcorr frame size:%d start lag:%d nlags:%d\n",
            session->dp.size_cir_buffer, session->dp.size, session->dp.start, session->dp.nlags);
  }

  session->dp.num_active_frames = 0;
  session->dp.first_time = 1;

  return 0;
}


/*--------------------------------------------------------------------*/
int dp_f0(
        get_f0_session *session,
        float* fdata,
        int buff_size,
        int sdstep,
        double freq,
        float **f0p_pt,
        float **vuvp_pt,
        float **rms_speech_pt,
        float **acpkp_pt,
        int *vecsize,
        int last_time
        )
{
    float  maxval, engref, *sta, *rms_ratio, *dsdata;
    register float ttemp, ftemp, ft1, ferr, err, errmin;
    register int  i, j, k, loc1, loc2;
    int   nframes, maxloc, ncand, ncandp, minloc,
            decimate, samsds;

    Stat *stat = NULL;
    F0_params *par = session->par;

    nframes = get_Nframes((long) buff_size, session->dp.pad, session->dp.step); /* # of whole frames */

    if(debug_level)
        Fprintf(stderr,
                "%s: ******* Computing %d dp frames ******** from %d points\n",
                ProgName, nframes, buff_size);

    /* Now downsample the signal for coarse peak estimates. */

    decimate = freq/2000.0;     /* downsample to about 2kHz */
    if (decimate <= 1)
        dsdata = fdata;
    else {
        samsds = ((nframes-1) * session->dp.step + session->dp.ncomp) / decimate;
        dsdata = downsample(
                    session,
                    fdata,
                    buff_size,
                    sdstep,
                    freq,
                    &samsds,
                    decimate,
                    last_time
                    );
        if (!dsdata) {
            Fprintf(stderr, "%s: can't get downsampled data.\n", ProgName);
            return 1;
        }
    }

    /* Get a function of the "stationarity" of the speech signal. */

    stat = get_stationarity(session, fdata, freq, buff_size, nframes, session->dp.step, session->dp.first_time);
    if (!stat) {
        Fprintf(stderr, "%s: can't get stationarity\n", ProgName);
        return(1);
    }
    sta = stat->stat;
    rms_ratio = stat->rms_ratio;

    /***********************************************************************/
    /* MAIN FUNDAMENTAL FREQUENCY ESTIMATION LOOP */
    /***********************************************************************/
    if(!session->dp.first_time && nframes > 0)
        session->dp.headF = session->dp.headF->next;

    for(i = 0; i < nframes; i++) {

        /* NOTE: This buffer growth provision is probably not necessary.
       It was put in (with errors) by Derek Lin and apparently never
       tested.  My tests and analysis suggest it is completely
       superfluous. DT 9/5/96 */
        /* Dynamically allocating more space for the circular buffer */
        if(session->dp.headF == session->dp.tailF->prev){
            Frame_getf0 *frm;

            if(session->dp.cir_buff_growth_count > 5){
                Fprintf(stderr,
                        "%s: too many requests (%d) for dynamically allocating space.\n   There may be a problem in finding converged path.\n",
                        ProgName, session->dp.cir_buff_growth_count);
                return(1);
            }
            if(debug_level)
                Fprintf(stderr, "%s: allocating %d more frames for DP circ. buffer.\n",
                        ProgName, session->dp.size_cir_buffer);
            frm = alloc_frame(session->dp.nlags, par->n_cands);
            session->dp.headF->next = frm;
            frm->prev = session->dp.headF;
            for(k=1; k<session->dp.size_cir_buffer; k++){
                frm->next = alloc_frame(session->dp.nlags, par->n_cands);
                frm->next->prev = frm;
                frm = frm->next;
            }
            frm->next = session->dp.tailF;
            session->dp.tailF->prev = frm;
            session->dp.cir_buff_growth_count++;
        }

        session->dp.headF->rms = stat->rms[i];
        if (get_fast_cands(session, fdata, dsdata, i,
                           session->dp.step, session->dp.size,
                           decimate, session->dp.start, session->dp.nlags,
                           &engref, &maxloc, &maxval,
                           session->dp.headF->cp, session->dp.peaks,
                           session->dp.locs, &ncand) != 0) {
            return 1; // failure
        }


        /*    Move the peak value and location arrays into the dp structure */
        {
            register float *ftp1, *ftp2;
            register short *sp1;
            register int *sp2;

            for(ftp1 = session->dp.headF->dp->pvals, ftp2 = session->dp.peaks,
                sp1 = session->dp.headF->dp->locs, sp2 = session->dp.locs, j=ncand; j--; ) {
                *ftp1++ = *ftp2++;
                *sp1++ = *sp2++;
            }
            *sp1 = -1;          /* distinguish the UNVOICED candidate */
            *ftp1 = maxval;
            session->dp.headF->dp->mpvals[ncand] = session->dp.vbias+maxval; /* (high cost if cor. is high)*/
        }

        /* Apply a lag-dependent weight to the peaks to encourage the selection
       of the first major peak.  Translate the modified peak values into
       costs (high peak ==> low cost). */
        for(j=0; j < ncand; j++){
            ftemp = 1.0 - ((float)session->dp.locs[j] * session->dp.lagwt);
            session->dp.headF->dp->mpvals[j] = 1.0 - (session->dp.peaks[j] * ftemp);
        }
        ncand++;                        /* include the unvoiced candidate */
        session->dp.headF->dp->ncands = ncand;

        /*********************************************************************/
        /*    COMPUTE THE DISTANCE MEASURES AND ACCUMULATE THE COSTS.       */
        /*********************************************************************/

        ncandp = session->dp.headF->prev->dp->ncands;
        for(k=0; k<ncand; k++){ /* for each of the current candidates... */
            minloc = 0;
            errmin = FLT_MAX;
            if((loc2 = session->dp.headF->dp->locs[k]) > 0) { /* current cand. is voiced */
                for(j=0; j<ncandp; j++){ /* for each PREVIOUS candidate... */
                    /*    Get cost due to inter-frame period change. */
                    loc1 = session->dp.headF->prev->dp->locs[j];
                    if (loc1 > 0) { /* prev. was voiced */
                        ftemp = log(((double) loc2) / loc1);
                        ttemp = fabs(ftemp);
                        ft1 = session->dp.fdouble + fabs(ftemp + session->dp.ln2);
                        if (ttemp > ft1)
                            ttemp = ft1;
                        ft1 = session->dp.fdouble + fabs(ftemp - session->dp.ln2);
                        if (ttemp > ft1)
                            ttemp = ft1;
                        ferr = ttemp * session->dp.freqwt;
                    } else {            /* prev. was unvoiced */
                        ferr = session->dp.tcost + (session->dp.tfact_s * sta[i]) + (session->dp.tfact_a / rms_ratio[i]);
                    }
                    /*    Add in cumulative cost associated with previous peak. */
                    err = ferr + session->dp.headF->prev->dp->dpvals[j];
                    if(err < errmin){   /* find min. cost */
                        errmin = err;
                        minloc = j;
                    }
                }
            } else {                    /* this is the unvoiced candidate */
                for(j=0; j<ncandp; j++){ /* for each PREVIOUS candidate... */

                    /*    Get voicing transition cost. */
                    if (session->dp.headF->prev->dp->locs[j] > 0) { /* previous was voiced */
                        ferr = session->dp.tcost + (session->dp.tfact_s * sta[i]) + (session->dp.tfact_a * rms_ratio[i]);
                    }
                    else
                        ferr = 0.0;
                    /*    Add in cumulative cost associated with previous peak. */
                    err = ferr + session->dp.headF->prev->dp->dpvals[j];
                    if(err < errmin){   /* find min. cost */
                        errmin = err;
                        minloc = j;
                    }
                }
            }
            /* Now have found the best path from this cand. to prev. frame */
            if (session->dp.first_time && i==0) {               /* this is the first frame */
                session->dp.headF->dp->dpvals[k] = session->dp.headF->dp->mpvals[k];
                session->dp.headF->dp->prept[k] = 0;
            } else {
                session->dp.headF->dp->dpvals[k] = errmin + session->dp.headF->dp->mpvals[k];
                session->dp.headF->dp->prept[k] = minloc;
            }
        } /*    END OF THIS DP FRAME */

        if (i < nframes - 1)
            session->dp.headF = session->dp.headF->next;

        if (debug_level >= 2) {
            Fprintf(stderr,"%d engref:%10.0f max:%7.5f loc:%4d\n",
                    i,engref,maxval,maxloc);
        }

    } /* end for (i ...) */

    /***************************************************************/
    /* DONE WITH FILLING DP STRUCTURES FOR THE SET OF SAMPLED DATA */
    /*    NOW FIND A CONVERGED DP PATH                             */
    /***************************************************************/

    *vecsize = 0;                       /* # of output frames returned */

    session->dp.num_active_frames += nframes;

    if( session->dp.num_active_frames >= session->dp.size_frame_hist  || last_time ){
        Frame_getf0 *frm;
        int  num_paths, best_cand, frmcnt, checkpath_done = 1;
        float patherrmin;

        if(debug_level)
            Fprintf(stderr, "%s: available frames for backtracking: %d\n",
                    ProgName, session->dp.num_active_frames);

        patherrmin = FLT_MAX;
        best_cand = 0;
        num_paths = session->dp.headF->dp->ncands;

        /* Get the best candidate for the final frame and initialize the
       paths' backpointers. */
        frm = session->dp.headF;
        for(k=0; k < num_paths; k++) {
            if (patherrmin > session->dp.headF->dp->dpvals[k]){
                patherrmin = session->dp.headF->dp->dpvals[k];
                best_cand = k;  /* index indicating the best candidate at a path */
            }
            session->dp.pcands[k] = frm->dp->prept[k];
        }

        if(last_time){     /* Input data was exhausted. force final outputs. */
            session->dp.cmpthF = session->dp.headF;             /* Use the current frame as starting point. */
        } else {
            /* Starting from the most recent frame, trace back each candidate's
  best path until reaching a common candidate at some past frame. */
            frmcnt = 0;
            while (1) {
                frm = frm->prev;
                frmcnt++;
                checkpath_done = 1;
                for(k=1; k < num_paths; k++){ /* Check for convergence. */
                    if(session->dp.pcands[0] != session->dp.pcands[k])
                        checkpath_done = 0;
                }
                if( ! checkpath_done) { /* Prepare for checking at prev. frame. */
                    for(k=0; k < num_paths; k++){
                        session->dp.pcands[k] = frm->dp->prept[session->dp.pcands[k]];
                    }
                } else {        /* All paths have converged. */
                    session->dp.cmpthF = frm;
                    best_cand = session->dp.pcands[0];
                    if(debug_level)
                        Fprintf(stderr,
                                "%s: paths went back %d frames before converging\n",
                                ProgName, frmcnt);
                    break;
                }
                if(frm == session->dp.tailF){   /* Used all available data? */
                    if( session->dp.num_active_frames < session->dp.size_frame_out) { /* Delay some more? */
                        checkpath_done = 0; /* Yes, don't backtrack at this time. */
                        session->dp.cmpthF = NULL;
                    } else {            /* No more delay! Force best-guess output. */
                        checkpath_done = 1;
                        session->dp.cmpthF = session->dp.headF;
                        Fprintf(stderr,
                                "%s: WARNING: no converging path found after going back %d frames, will use the lowest cost path\n",
                                ProgName, session->dp.num_active_frames);
                    }
                    break;
                } /* end if (frm ...) */
            }   /* end while (1) */
        } /* end if (last_time) ... else */

        /*************************************************************/
        /* BACKTRACKING FROM cmpthF (best_cand) ALL THE WAY TO tailF    */
        /*************************************************************/
        i = 0;
        frm = session->dp.cmpthF;       /* Start where convergence was found (or faked). */
        while( frm != session->dp.tailF->prev && checkpath_done){
            if( i == session->dp.output_buf_size ){ /* Need more room for outputs? */
                session->dp.output_buf_size *= 2;
                if(debug_level)
                    Fprintf(stderr,
                            "%s: reallocating space for output frames: %d\n",
                            ProgName, session->dp.output_buf_size);
                session->dp.rms_speech = (float *) realloc((char *) session->dp.rms_speech,
                                                           sizeof(float) * session->dp.output_buf_size);
                if (!session->dp.rms_speech) {
                    fprintf(stderr, "rms_speech realloc failed in dp_f0()\n");
                    return 1;
                }
                session->dp.f0p = (float *) realloc((char *) session->dp.f0p,
                                                    sizeof(float) * session->dp.output_buf_size);
                if (!session->dp.f0p) {
                    fprintf(stderr, "f0p realloc failed in dp_f0()\n");
                    return 1;
                }
                session->dp.vuvp = (float *) realloc(session->dp.vuvp, sizeof(float) * session->dp.output_buf_size);
                if (!session->dp.vuvp) {
                    fprintf(stderr, "vuvp realloc failed in dp_f0()\n");
                    return 1;
                }
                session->dp.acpkp = (float *) realloc(session->dp.acpkp, sizeof(float) * session->dp.output_buf_size);
                if (!session->dp.acpkp) {
                    fprintf(stderr, "acpkp realloc failed in dp_f0()\n");
                    return 1;
                }
            }
            session->dp.rms_speech[i] = frm->rms;
            session->dp.acpkp[i] =  frm->dp->pvals[best_cand];
            loc1 = frm->dp->locs[best_cand];
            session->dp.vuvp[i] = 1.0;
            best_cand = frm->dp->prept[best_cand];
            ftemp = loc1;
            if(loc1 > 0) {              /* Was f0 actually estimated for this frame? */
                if (loc1 > session->dp.start && loc1 < session->dp.stop) { /* loc1 must be a local maximum. */
                    float cormax, cprev, cnext, den;

                    j = loc1 - session->dp.start;
                    cormax = frm->cp->correl[j];
                    cprev = frm->cp->correl[j+1];
                    cnext = frm->cp->correl[j-1];
                    den = (2.0 * ( cprev + cnext - (2.0 * cormax) ));
                    /*
    * Only parabolic interpolate if cormax is indeed a local
    * turning point. Find peak of curve that goes though the 3 points
    */

                    if (fabs(den) > 0.000001)
                        ftemp += 2.0 - ((((5.0*cprev)+(3.0*cnext)-(8.0*cormax))/den));
                }
                session->dp.f0p[i] = freq/ftemp;
            } else {            /* No valid estimate; just fake some arbitrary F0. */
                session->dp.f0p[i] = 0;
                session->dp.vuvp[i] = 0.0;
            }
            frm = frm->prev;

            if (debug_level >= 2)
                Fprintf(stderr," i:%4d%8.1f%8.1f\n",i,session->dp.f0p[i],session->dp.vuvp[i]);
            /* f0p[i] starts from the most recent one */
            /* Need to reverse the order in the calling function */
            i++;
        } /* end while() */
        if (checkpath_done){
            *vecsize = i;
            session->dp.tailF = session->dp.cmpthF->next;
            session->dp.num_active_frames -= *vecsize;
        }
    } /* end if() */

    if (debug_level)
        Fprintf(stderr, "%s: writing out %d frames.\n", ProgName, *vecsize);

    *f0p_pt = session->dp.f0p;
    *vuvp_pt = session->dp.vuvp;
    *acpkp_pt = session->dp.acpkp;
    *rms_speech_pt = session->dp.rms_speech;

    if(session->dp.first_time) session->dp.first_time = 0;
    return(0);
}


/*--------------------------------------------------------------------*/
static void free_frame(Frame_getf0 *frm)
{
    if (frm) {
        if (frm->dp) {
            if (frm->dp->locs)
                free(frm->dp->locs);
            if (frm->dp->pvals)
                free(frm->dp->pvals);
            if (frm->dp->mpvals)
                free(frm->dp->mpvals);
            if (frm->dp->prept)
                free(frm->dp->prept);
            if (frm->dp->dpvals)
                free(frm->dp->dpvals);
            free(frm->dp);
        }
        if (frm->cp) {
            if (frm->cp->correl)
                free(frm->cp->correl);
            free(frm->cp);
        }
        free(frm);
    }
}

Frame_getf0 *alloc_frame(int nlags, int ncands)
{
    Frame_getf0 *frm;
    int j;

    frm = (Frame_getf0 *) calloc(1, sizeof(Frame_getf0));
    if (!frm) {
        fprintf(stderr, "malloc failed in alloc_frame\n");
        return 0;
    }
    frm->dp = (Dprec *) calloc(1, sizeof(Dprec));
    if (!frm->dp) {
        fprintf(stderr, "frm->dp malloc failed in alloc_frame\n");
        free_frame(frm);
        return 0;
    }
    frm->dp->ncands = 0;
    frm->cp = (Cross *) calloc(1, sizeof(Cross));
    if (!frm->cp) {
        fprintf(stderr, "frm->cp malloc failed in alloc_frame\n");
        free_frame(frm);
        return 0;
    }
    frm->cp->correl = (float *) malloc(sizeof(float) * nlags);
    if (!frm->cp->correl) {
        fprintf(stderr, "frm->cp->correl malloc failed\n");
        free_frame(frm);
        return 0;
    }
    /* Allocate space for candidates and working arrays. */
    frm->dp->locs = (short*)malloc(sizeof(short) * ncands);
    if (!frm->dp->locs) {
        fprintf(stderr, "frm->dp->locs malloc failed in alloc_frame()\n");
        free_frame(frm);
        return 0;
    }
    frm->dp->pvals = (float*)malloc(sizeof(float) * ncands);
    if (!frm->dp->pvals) {
        fprintf(stderr, "frm->dp->pvals malloc failed in alloc_frame()\n");
        free_frame(frm);
        return 0;
    }
    frm->dp->mpvals = (float*)malloc(sizeof(float) * ncands);
    if (!frm->dp->mpvals) {
        fprintf(stderr, "frm->dp->mpvals malloc failed in alloc_frame()\n");
        free_frame(frm);
        return 0;
    }
    frm->dp->prept = (short*)malloc(sizeof(short) * ncands);
    if (!frm->dp->prept) {
        fprintf(stderr, "frm->dp->prept malloc failed in alloc_frame()\n");
        free_frame(frm);
        return 0;
    }
    frm->dp->dpvals = (float*)malloc(sizeof(float) * ncands);
    if (!frm->dp->dpvals) {
        fprintf(stderr, "frm->dp->dpvals malloc failed in alloc_frame()\n");
        free_frame(frm);
        return 0;
    }

    /*  Initialize the cumulative DP costs to zero */
    for(j = ncands-1; j >= 0; j--)
        frm->dp->dpvals[j] = 0.0;

    return(frm);
}


/*--------------------------------------------------------------------*/
/* push window stat to stack, and pop the oldest one */

static int save_windstat(
              get_f0_session *session,
              float *rho,
              int order,
              float err,
              float rms
              )
{
    int i,j;

    if(session->dp.wReuse > 1){               /* push down the stack */
        for(j=1; j<session->dp.wReuse; j++){
            for(i=0;i<=order; i++) session->dp.windstat[j-1].rho[i] = session->dp.windstat[j].rho[i];
            session->dp.windstat[j-1].err = session->dp.windstat[j].err;
            session->dp.windstat[j-1].rms = session->dp.windstat[j].rms;
        }
        for(i=0;i<=order; i++) session->dp.windstat[session->dp.wReuse-1].rho[i] = rho[i]; /*save*/
        session->dp.windstat[session->dp.wReuse-1].err = err;
        session->dp.windstat[session->dp.wReuse-1].rms = rms;
        return 1;
    } else if (session->dp.wReuse == 1) {
        for(i=0;i<=order; i++) session->dp.windstat[0].rho[i] = rho[i];  /* save */
        session->dp.windstat[0].err = err;
        session->dp.windstat[0].rms = rms;
        return 1;
    } else
        return 0;
}


/*--------------------------------------------------------------------*/
static int retrieve_windstat(
        get_f0_session *session,
        float *rho,
        int order,
        float *err,
        float *rms
        )
{
    Windstat wstat;
    int i;
       
    if(session->dp.wReuse){
        wstat = session->dp.windstat[0];
        for(i=0; i<=order; i++) rho[i] = wstat.rho[i];
        *err = wstat.err;
        *rms = wstat.rms;
        return 1;
    }
    else return 0;
}


/*--------------------------------------------------------------------*/
static float get_similarity(
        get_f0_session *session,
        int order,
        int size,
        float *pdata,
        float *cdata,
        float *rmsa,
        float *rms_ratio,
        float pre,
        float stab,
        int w_type,
        int init
        )
{
    float rho3[BIGSORD+1], err3, rms3, rmsd3, b0, t, a2[BIGSORD+1],
            rho1[BIGSORD+1], a1[BIGSORD+1], b[BIGSORD+1], err1, rms1, rmsd1;

    /* (In the lpc() calls below, size-1 is used, since the windowing and
   preemphasis function assumes an extra point is available in the
   input data array.  This condition is apparently no longer met after
   Derek's modifications.) */

    /* get current window stat */
    lpc(session, order, stab, size-1, cdata,
        a2, rho3, (float *) NULL, &err3, &rmsd3, pre, w_type);
    rms3 = wind_energy(session, cdata, size, w_type);

    if(!init) {
        /* get previous window stat */
        if( !retrieve_windstat(session, rho1, order, &err1, &rms1)){
            lpc(session, order, stab, size-1, pdata,
                a1, rho1, (float *) NULL, &err1, &rmsd1, pre, w_type);
            rms1 = wind_energy(session, pdata, size, w_type);
        }
        a_to_aca(a2+1,b,&b0,order);
        t = itakura(order,b,&b0,rho1+1,&err1) - .8;
        if(rms1 > 0.0)
            *rms_ratio = (0.001 + rms3)/rms1;
        else
            if(rms3 > 0.0)
                *rms_ratio = 2.0;       /* indicate some energy increase */
            else
                *rms_ratio = 1.0;       /* no change */
    } else {
        *rms_ratio = 1.0;
        t = 10.0;
    }
    *rmsa = rms3;
    save_windstat(session, rho3, order, err3, rms3);
    return((float)(0.2/t));
}


/* -------------------------------------------------------------------- */
/* This is an ad hoc signal stationarity function based on Itakura
 * distance and relative amplitudes.
 */
/*
  This illustrates the window locations when the very first frame is read.
  It shows an example where each frame step |  .  | is 10 msec.  The
  frame step size is variable.  The window size is always 30 msec.
  The window centers '*' is always 20 msec apart.
  The windows cross each other right at the center of the DP frame, or
  where the '.' is.

                          ---------*---------   current window

              ---------*---------  previous window

  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |
              ^           ^  ^
              ^           ^  ^
              ^           ^  fdata
              ^           ^
              ^           q
              p

                          ---
                          ind

  fdata, q, p, ind, are variables used below.

*/

static void free_stat(Stat* stat)
{
    if (stat) {
        if (stat->stat)
            free((char *) stat->stat);
        if (stat->rms)
            free((char *) stat->rms);
        if (stat->rms_ratio)
            free((char *) stat->rms_ratio);
        free((char *) stat);
    }
}

static Stat* get_stationarity(
        get_f0_session *session,
        float *fdata,
        double freq,
        int buff_size,
        int nframes,
        int frame_step,
        int first_time
        )
{
    float preemp = 0.4, stab = 30.0;
    float *p, *q, *r, *datend;
    int ind, i, j, m, size, order, agap, w_type = 3;

    agap = (int) (STAT_AINT *freq);
    size = (int) (STAT_WSIZE * freq);
    ind = (agap - size) / 2;

    if( session->dp.nframes_old < nframes || !session->dp.stat || first_time){
        /* move this to init_dp_f0() later */
        session->dp.nframes_old = nframes;
        if (session->dp.stat)
            free_stat(session->dp.stat);
        session->dp.stat = (Stat *) calloc(1, nframes * sizeof(Stat));
        if (!session->dp.stat) {
            fprintf(stderr, "stat malloc failed in get_stationarity\n");
            return 0;
        }
        session->dp.stat->stat = (float*)malloc(sizeof(float)*nframes);
        if (!session->dp.stat->stat) {
            fprintf(stderr, "stat->stat malloc failed in get_stationarity\n");
            free_stat(session->dp.stat);
            return 0;
        }
        session->dp.stat->rms = (float*)malloc(sizeof(float)*nframes);
        if (!session->dp.stat->rms) {
            fprintf(stderr, "stat->rms malloc failed in get_stationarity\n");
            free_stat(session->dp.stat);
            return 0;
        }
        session->dp.stat->rms_ratio = (float*)malloc(sizeof(float)*nframes);
        if (!session->dp.stat->rms_ratio) {
            fprintf(stderr, "stat->ratio malloc failed in get_stationarity\n");
            free_stat(session->dp.stat);
            return 0;
        }
        session->dp.memsize = (int) (STAT_WSIZE * freq) + (int) (STAT_AINT * freq);
        session->dp.mem = (float *) malloc( sizeof(float) * session->dp.memsize);
        if (!session->dp.mem) {
            fprintf(stderr, "mem malloc failed in get_stationarity()\n");
            return 0;
        }
        for(j=0; j<session->dp.memsize; j++) session->dp.mem[j] = 0;
    }

    if(nframes == 0) return(session->dp.stat);

    q = fdata + ind;
    datend = fdata + buff_size;

    if((order = 2.0 + (freq/1000.0)) > BIGSORD) {
        Fprintf(stderr,
                "%s: Optimim order (%d) exceeds that allowable (%d); reduce Fs\n",
                ProgName, order, BIGSORD);
        order = BIGSORD;
    }

    /* prepare for the first frame */
    for(j=session->dp.memsize/2, i=0; j<session->dp.memsize; j++, i++) session->dp.mem[j] = fdata[i];

    /* never run over end of frame, should already taken care of when read */

    for(j=0, p = q - agap; j < nframes; j++, p += frame_step, q += frame_step){
        if( (p >= fdata) && (q >= fdata) && ( q + size <= datend) )
            session->dp.stat->stat[j] = get_similarity(
                        session,
                        order,
                        size,
                        p,
                        q,
                        &(session->dp.stat->rms[j]),
                        &(session->dp.stat->rms_ratio[j]),
                        preemp,
                        stab,
                        w_type,
                        0);
        else {
            if(first_time) {
                if( (p < fdata) && (q >= fdata) && (q+size <=datend) )
                    session->dp.stat->stat[j] = get_similarity(
                                session,
                                order,
                                size,
                                NULL,
                                q,
                                &(session->dp.stat->rms[j]),
                                &(session->dp.stat->rms_ratio[j]),
                                preemp,
                                stab,
                                w_type,
                                1);
                else{
                    session->dp.stat->rms[j] = 0.0;
                    session->dp.stat->stat[j] = 0.01 * 0.2;   /* a big transition */
                    session->dp.stat->rms_ratio[j] = 1.0;   /* no amplitude change */
                }
            } else {
                if( (p<fdata) && (q+size <=datend) ){
                    session->dp.stat->stat[j] = get_similarity(
                                session,
                                order,
                                size,
                                session->dp.mem,
                                session->dp.mem + (session->dp.memsize/2) + ind,
                                &(session->dp.stat->rms[j]),
                                &(session->dp.stat->rms_ratio[j]),
                                preemp,
                                stab,
                                w_type,
                                0);
                    /* prepare for the next frame_step if needed */
                    if(p + frame_step < fdata ){
                        for( m=0; m<(session->dp.memsize-frame_step); m++)
                            session->dp.mem[m] = session->dp.mem[m+frame_step];
                        r = q + size;
                        for( m=0; m<frame_step; m++)
                            session->dp.mem[session->dp.memsize-frame_step+m] = *r++;
                    }
                }
            }
        }
    }

    /* last frame, prepare for next call */
    for(j=(session->dp.memsize/2)-1, p=fdata + (nframes * frame_step)-1; j>=0; j-- )
        session->dp.mem[j] = *p--;
    return(session->dp.stat);
}


/* -------------------------------------------------------------------- */
/*      Round the argument to the nearest integer.                      */

int Round(double flnum)
{
    return((flnum >= 0.0) ? (int)(flnum + 0.5) : (int)(flnum - 0.5));
}


/// Free up dynamic memories allocated by dp_f0 module
void close_dp_f0(get_f0_session *session)
{
    if (session->dp.headF) {
        Frame_getf0 *p = session->dp.headF;
        Frame_getf0 *q;
        while (p) {
            q = p->next;
            free_frame(p);
            p = q;
            if (p == session->dp.headF)
                break;
        }
    }
    if (session->dp.pcands)
        free(session->dp.pcands);
    if (session->dp.rms_speech)
        free(session->dp.rms_speech);
    if (session->dp.f0p)
        free(session->dp.f0p);
    if (session->dp.vuvp)
        free(session->dp.vuvp);
    if (session->dp.acpkp)
        free(session->dp.acpkp);
    if (session->dp.peaks)
        free(session->dp.peaks);
    if (session->dp.locs)
        free(session->dp.locs);
    if (session->dp.wReuse && session->dp.windstat)
        free(session->dp.windstat);
    if (session->dp.stat)
        free_stat(session->dp.stat);
    if (session->dp.mem)
        free(session->dp.mem);

    // reset values for future use
    session->dp.headF = 0;
    session->dp.tailF = 0;
    session->dp.pcands = 0;
    session->dp.rms_speech = 0;
    session->dp.f0p = 0;
    session->dp.vuvp = 0;
    session->dp.acpkp = 0;
    session->dp.peaks = 0;
    session->dp.locs = 0;
    session->dp.windstat = 0;
    session->dp.wReuse = 0;
    session->dp.first_time = 1;
    session->dp.cir_buff_growth_count = 0;
    session->dp.stat = 0;
    session->dp.nframes_old = 0;

    // get_cands
    if (session->gcs.foutput)
        free(session->gcs.foutput);
    if (session->gcs.co)
        free(session->gcs.co);
    if (session->gcs.mem)
        free(session->gcs.mem);

    // reset values for future use
    session->gcs.foutput = 0;
    session->gcs.ncoeff = 127;
    session->gcs.ncoefft = 0;
    session->gcs.co = 0;
    session->gcs.mem = 0;
    session->gcs.fsize = 0;
    session->gcs.resid = 0;

    // sigproc
    if (session->sp.din)
        free(session->sp.din);
    if (session->sp.c_wind)
        free(session->sp.c_wind);
    if (session->sp.h_wind)
        free(session->sp.h_wind);
    if (session->sp.hn_wind)
        free(session->sp.hn_wind);
    if (session->sp.we_dwind)
        free(session->sp.we_dwind);
    if (session->sp.lpc_dwind)
        free(session->sp.lpc_dwind);
    if (session->sp.f_dbdata)
        free(session->sp.f_dbdata);
    if (session->sp.fi_dbdata)
        free(session->sp.fi_dbdata);
    session->sp.din = 0;
    session->sp.c_wind = 0;
    session->sp.h_wind = 0;
    session->sp.hn_wind = 0;
    session->sp.we_dwind = 0;
    session->sp.lpc_dwind = 0;
    session->sp.f_dbdata = 0;
    session->sp.fi_dbdata = 0;
    session->sp.n0 = 0;
    session->sp.c_wsize = 0;
    session->sp.h_wsize = 0;
    session->sp.hn_wsize = 0;
    session->sp.we_nwind = 0;
    session->sp.lpc_nwind = 0;
    session->sp.f_dbsize = 0;
    session->sp.fi_dbsize = 0;
}
