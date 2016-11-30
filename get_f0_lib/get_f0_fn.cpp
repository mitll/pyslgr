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
 * Written by:  Derek Lin
 * Checked by:
 * Revised by:  David Talkin
 *
 * Brief description:  Estimates F0 using normalized cross correlation and
 *   dynamic programming.
 *
 */

#include "stdafx.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "speech_tools.h"

static int get_samples (Signal &x, int start, vec<float> &out);

#include "f0.h"

/*
 * Some consistency checks on parameter values.
 * Return a positive integer if any errors detected, 0 if none.
 */

int check_f0_params(F0_params *par, double sample_freq)
{
  int     error = 0;
  double  dstep;

  if((par->cand_thresh < 0.01) || (par->cand_thresh > 0.99)) {
    fprintf(stderr,
            "ERROR: cand_thresh parameter must be between [0.01, 0.99].\n");
    error++;
  }
  if((par->wind_dur > .1) || (par->wind_dur < .0001)) {
    fprintf(stderr,
            "ERROR: wind_dur parameter must be between [0.0001, 0.1].\n");
    error++;
  }
  if((par->n_cands > 100) || (par->n_cands < 3)){
    fprintf(stderr,
            "ERROR: n_cands parameter must be between [3,100].\n");
    error++;
  }
  if((par->max_f0 <= par->min_f0) || (par->max_f0 >= (sample_freq/2.0)) ||
     (par->min_f0 < (sample_freq/10000.0))){
    fprintf(stderr,
            "ERROR: min(max)_f0 parameter inconsistent with sampling frequency.\n");
    error++;
  }
  dstep = ((double)((int)(0.5 + (sample_freq * par->frame_step))))/sample_freq;
  if(dstep != par->frame_step) {
    fprintf(stderr,
            "Frame step set to %f to exactly match signal sample rate.\n",
            dstep);
    par->frame_step = dstep;
  }
  if((par->frame_step > 0.1) || (par->frame_step < (1.0/sample_freq))){
    fprintf(stderr,
            "ERROR: frame_step parameter must be between [1/sampling rate, 0.1].\n");
    error++;
  }

  return(error);
}

F0_params* get_f0_params()
{
    F0_params *par = (F0_params *) malloc(sizeof(F0_params));
    par->cand_thresh = 0.3;
    par->lag_weight = 0.3;
    par->freq_weight = 0.02;
    par->trans_cost = 0.005;
    par->trans_amp = 0.5;
    par->trans_spec = 0.5;
    par->voice_bias = 0.0;
    par->double_cost = 0.35;
    par->min_f0 = 50.0;
    par->max_f0 = 550.0;
    par->frame_step = 0.01;
    par->wind_dur = 0.0075;
    par->n_cands = 20;
    par->mean_f0 = 200.0;     /* unused */
    par->mean_f0_weight = 0.0;  /* unused */
    par->conditioning = 0;    /*unused */
    return par;
}


extern int init_dp_f0(get_f0_session*, double, long*, long*);
extern void close_dp_f0(get_f0_session*);
extern int dp_f0(get_f0_session*, float*, int, int, double, float**, float**, float**, float**, int*, int);


get_f0_session* init_get_f0()
{
    get_f0_session *session = (get_f0_session*) calloc(1, sizeof(get_f0_session));
    session->par = get_f0_params();
    session->dp.first_time = 1;
    session->gcs.ncoeff = 127;
    return session;
}

void close_get_f0(get_f0_session *session)
{
  free(session->par);
  free(session);
}

int get_f0_esps (Signal &x, get_f0_session *session, int start_sample, int end_sample, vector<float> &f0_store)
{
    long total_samps;
    long min_samps;
    long buff_size;
    long sdstep;
    long actsize;
    int done;
    float *f0p;
    float *vuvp;
    float *rms_speech;
    float *acpkp;
    int vecsize;

    if (end_sample < start_sample || start_sample < 0)
        return 1;

    double samplerate = x.sampling_freq();
    if (start_sample >= x.len())
        return 1;

    if (end_sample >= x.len())
        end_sample = x.len()-1;

    if (check_f0_params(session->par, samplerate))
        return 1;

    total_samps = end_sample - start_sample + 1;
    min_samps = (session->par->frame_step * 2.0 + session->par->wind_dur)*samplerate;
    if (total_samps < min_samps)
        return 1;

    if (init_dp_f0(session, samplerate, &buff_size, &sdstep) != 0)
        return 1;

    if (buff_size > INT_MAX || sdstep > INT_MAX) {
        close_dp_f0(session);
        return 1;
    }

    vec<float> fdata(buff_size);

    if (buff_size > total_samps)
        buff_size = total_samps;

    int curr_sample = start_sample, i;

    actsize = get_samples(x, curr_sample, fdata);

    while (1) {
        done = (actsize < buff_size) || (total_samps == buff_size);

        if (dp_f0(session, fdata.data, (int) actsize, (int) sdstep, samplerate, &f0p, &vuvp, &rms_speech, &acpkp, &vecsize, done) != 0) {
            close_dp_f0(session);
            return 1;
        }

        for (i=vecsize-1; i>=0; i--)
           f0_store.push_back(f0p[i]);

        if (done)
            break;

        curr_sample += sdstep;
        actsize = get_samples (x, curr_sample, fdata);

        total_samps -= sdstep;

        if (actsize > total_samps)
            actsize = total_samps;
    }

    close_dp_f0(session);

    return 0;
}

static int get_samples (Signal &x, int start, vec<float> &out)
{
   int buff_size, i, j, len;

   len = x.len();
   buff_size = out.len;

   for (i=0; i<buff_size; i++)
      out.data[i] = 0;

   for (i=start, j=0; (j<buff_size) && (i<len); i++, j++)
      out.data[j] = ((float) x.value[i]);

   return j;

}
