#ifndef F0_H
#define F0_H

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
 * @(#)f0.h     1.4 9/9/96 ERL
 * Brief description:
 *
 */

/* f0.h */
/* Some definitions used by the "Pitch Tracker Software". */

#include "f0_structs.h"
       
typedef struct f0_params {
float cand_thresh,      /* only correlation peaks above this are considered */
      lag_weight,       /* degree to which shorter lags are weighted */
      freq_weight,      /* weighting given to F0 trajectory smoothness */
      trans_cost,       /* fixed cost for a voicing-state transition */
      trans_amp,        /* amplitude-change-modulated VUV trans. cost */
      trans_spec,       /* spectral-change-modulated VUV trans. cost */
      voice_bias,       /* fixed bias towards the voiced hypothesis */
      double_cost,      /* cost for octave F0 jumps */
      mean_f0,          /* talker-specific mean F0 (Hz) */
      mean_f0_weight,   /* weight to be given to deviations from mean F0 */
      min_f0,           /* min. F0 to search for (Hz) */
      max_f0,           /* max. F0 to search for (Hz) */
      frame_step,       /* inter-frame-interval (sec) */
      wind_dur;         /* duration of correlation window (sec) */
int   n_cands,          /* max. # of F0 cands. to consider at each frame */
      conditioning;     /* Specify optional signal pre-conditioning. */
} F0_params;

/* Possible values returned by the function f0(). */
#define F0_OK           0
#define F0_NO_RETURNS   1
#define F0_TOO_FEW_SAMPLES      2
#define F0_NO_INPUT     3
#define F0_NO_PAR       4
#define F0_BAD_PAR      5
#define F0_BAD_INPUT    6
#define F0_INTERNAL_ERR 7

/* Bits to specify optional pre-conditioning of speech signals by f0() */
/* These may be OR'ed together to specify all preprocessing. */
#define F0_PC_NONE      0x00            /* no pre-processing */
#define F0_PC_DC        0x01            /* remove DC */
#define F0_PC_LP2000    0x02            /* 2000 Hz lowpass */
#define F0_PC_HP100     0x04            /* 100 Hz highpass */
#define F0_PC_AR        0x08            /* inf_order-order LPC inverse filter */
#define F0_PC_DIFF      0x010           /* 1st-order difference */


#define Fprintf (void)fprintf

F0_params* get_f0_params();

typedef struct {
    F0_params *par;

    // dp_f0
    struct {
        Frame_getf0 *headF;  /* current frame in the circular buffer */
        Frame_getf0 *tailF;  /* the frame where tracks start */
        Frame_getf0 *cmpthF; /* starting frame of converged path to backtrack */
        int *pcands;           /* array for backtracking in convergence check */
        int cir_buff_growth_count;
        int size_cir_buffer;   /* # of frames in circular DP buffer */
        int size_frame_hist;   /* # of frames required before convergence test */
        int size_frame_out;    /* # of frames before forcing output */
        int num_active_frames; /* # of frames from tailF to headF */
        int output_buf_size;   /* # of frames allocated to output buffers */

        /* DP parameters */
        float tcost;
        float tfact_a;
        float tfact_s;
        float frame_int;
        float vbias;
        float fdouble;
        float ln2;
        float freqwt;
        float lagwt;
        int step;
        int size;
        int nlags;
        int start;
        int stop;
        int ncomp;
        int *locs;
        short maxpeaks;
        int wReuse;            /* number of windows seen before resued */
        Windstat *windstat;
        float *f0p;
        float *vuvp;
        float *rms_speech;
        float *acpkp;
        float *peaks;
        int first_time;
        int pad;

        Stat *stat;
        int nframes_old;
        int memsize;
        float *mem;
    } dp;

    // get_cands
    struct {
        // downsample()
        float b[2048];
        float *foutput;
        int ncoeff;
        int ncoefft;

        // do_ffir()
        float *co;
        float *mem;
        float state[1000];
        int fsize;
        int resid;
    } gcs;

    // sigproc
    struct {
        // get_window()
        float *din;
        int n0;

        // cwindow()
        int c_wsize;
        float *c_wind;

        // hwindow()
        int h_wsize;
        float *h_wind;

        // hnwindow()
        int hn_wsize;
        float *hn_wind;

        // wind_energy()
        int we_nwind;
        float *we_dwind;

        // lpc()
        int lpc_nwind;
        float *lpc_dwind;

        // crossf()
        float *f_dbdata;
        int f_dbsize;

        // crossfi()
        float *fi_dbdata;
        int fi_dbsize;
    } sp;

} get_f0_session;

extern F0_params *new_f0_params();
extern int get_window();
extern int lpc(
    get_f0_session *session,
    int lpc_ord,
    float lpc_stabl,
    int wsize,
    float *data,
    float *lpca,
    float *ar,
    float *lpck,
    float *normerr,
    float *rms,
    float preemp,
    int type);
extern int Round(double flnum);
extern int window(
    get_f0_session *session,
    register float *din,
    register float *dout,
    register int n,
    register float preemp,
    int type);
extern void a_to_aca (float *a, float *b, float *c, register int p);
extern void cross(), autoc(), durbin();
extern int crossf(
        get_f0_session *session,
        float *data,
        int size,
        int start,
        int nlags,
        float *engref,
        int *maxloc,
        float *maxval,
        float *correl
        );
extern int crossfi(
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
        );


typedef void (*get_f0_callback_t)(float *f0p, float *vuvp, float *rms_speech, float *acpkp, int vecsize, void *user_data);

int get_f0_esps (Signal &x, get_f0_session *session, int start_sample, int end_sample, vector<float> &f0_store);
get_f0_session *init_get_f0();
void close_get_f0(get_f0_session *session);

#endif // F0_H
