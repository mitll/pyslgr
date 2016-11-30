//
// Copyright 2016 MIT Lincoln Laboratory, Massachusetts Institute of Technology
//
// Licensed under the Apache License, Version 2.0 (the "License"); you may not use these files except in compliance with
// the License.
//
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
// an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
// specific language governing permissions and limitations under the License.
//

//
// GMMSAD implementation for Signal class
//

//
// Standard headers
//
#include "stdafx.h"
#include "speech_tools.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
using boost::property_tree::ptree;
list<Time_Interval> get_tgt_segments (list<Time_Interval> &seg, int tgt_lbl);


//
// Local functions
//
static void find_intervals (vec<float> &frame_scores, int num_models, list<Time_Interval> &segments, float win_inc);
static void smooth_frame_scores (vec<float> &frame_scores, int num_models, int win_len);
static list<Time_Interval> smooth_segments(list<Time_Interval> &segments, float min_gap_dur, float min_seg_dur);

list<Time_Interval> Signal::gmmsad (Signal &x, string &feat_config, vector<GMM_model> &gmm_models, int label_window, int gmmsad_keep,
   float min_gap_dur, float min_seg_dur)
{
   // Extract features
	MFCC_Features f;
	f.process(x, feat_config);
	f.rasta();
   f.delta(2);
   f.set_outfeat("fd");
   float win_inc = 0.001*f.get_win_inc_ms();

   // Set up for scoring
   int num_models = (int) gmm_models.size();
   vec<float> frame_scores(num_models*f.num_vec());
   vec<float> frame_scores_one_model(f.num_vec());

   // Score GMM models independently
   int i, j;
   vector<GMM_model>::iterator it;
   float model_score; 
	vec<short> empty_frame_scores(0);

   for (it=gmm_models.begin(), i=0; it<gmm_models.end(); it++, i++) {
      model_score = gmm_models[i].score(f, frame_scores_one_model, empty_frame_scores);
      for (j=0; j<f.num_vec(); j++) 
         frame_scores.data[i+j*num_models] = frame_scores_one_model.data[j];
   }

   // Find segments
   list<Time_Interval> segments;
   smooth_frame_scores(frame_scores, num_models, label_window);
   find_intervals(frame_scores, num_models, segments, win_inc);

   // Get target segments
   list<Time_Interval> segments_tgt = get_tgt_segments(segments, gmmsad_keep);

   // Smooth segments
   list<Time_Interval> segments_smoothed = smooth_segments(segments_tgt, min_gap_dur, min_seg_dur);

   // Final result
   return segments_smoothed;
     
} // gmmsad

Time_Interval::Time_Interval (float start_in, float end_in, int lbl_in, float score_in)
{
   start = start_in;
   end = end_in;
   label = lbl_in;
   score = score_in;
}

list<Time_Interval> get_tgt_segments (list<Time_Interval> &seg, int tgt_lbl)
{
   list<Time_Interval> out;
   list<Time_Interval>::iterator it;

   for (it=seg.begin(); it!=seg.end(); it++) {
      if (it->label==tgt_lbl)
         out.push_back(*it);
   }
   return out;

}

list<Time_Interval> combine_segments (list<Time_Interval> &s1, list<Time_Interval> &s2)
{
   // compute s1 and not s2
   // if s2 is empty, nothing to remove from s1
   if (s2.size()==0)
      return s1;

   list<Time_Interval> output = s1;
   list<Time_Interval>::iterator it1, it2;
   it1 = output.begin();
   it2 = s2.begin();

   float t1_start, t1_end;
   float t2_start, t2_end;
   int tmp_lbl;

   while (it1!=output.end() && it2!=s2.end()) {
      t1_start = it1->start;
      t1_end = it1->end;
      t2_start = it2->start;
      t2_end = it2->end;

      if (t1_start >= t2_end) {
         it2++;
      } else if (t1_end <= t2_start) {
         it1++;
      } else if ((t1_start<=t2_start) && (t1_end>=t2_end)) { 
         // case 1 -- s1 completely contains s2
         tmp_lbl = it1->label;
         it1->end = it2->start;
         it1++;
         output.insert(it1, Time_Interval(it2->end, t1_end, tmp_lbl));
         it1--;
         it2++;
      } else if ((t1_start<=t2_start) && (t1_end < t2_end)) {
         // case 2 -- right side of s1 overlaps with left part of s2
         it1->end = t2_start;
         it1++;
      } else if ((t1_start>t2_start) && (t1_end>=t2_end)) {
         // case 3 -- left side of s1 overlaps with right part of s2
         it1->start = t2_end;
         it2++;
      } else if ((t1_start>=t2_start) && (t1_end<=t2_end)) {
         it1 = output.erase(it1);
      } else {
         throw ST_exception("Combine marks failed.");
      }
   }

   return output;

}

static list<Time_Interval> smooth_segments(list<Time_Interval> &segments, float min_gap_dur, float min_seg_dur)
{
   list<Time_Interval> seg_out;
   list<Time_Interval>::iterator it;

   float t1, t2;
   int lbl;
   for (it=segments.begin() ;it!=segments.end(); it++) {
      if (it==segments.begin()) {
         t1 = it->start;
         t2 = it->end;
         lbl = it->label;
      } else {
         if ((it->start-t2) < min_gap_dur) { // greedily combine segments
            t2 = it->end;
            continue;
         }
         if ((t2-t1) >= min_seg_dur)
            seg_out.push_back(Time_Interval(t1, t2, lbl));
         t1 = it->start;
         t2 = it->end;
      }
   }

   if ((t2-t1) > min_seg_dur)
      seg_out.push_back(Time_Interval(t1, t2, lbl));

   return seg_out;
}

static void find_intervals (vec<float> &frame_scores, int num_models, list<Time_Interval> &segments, float win_inc)
{
   int i, i_start, j, i_mx, curr_label;
   int num_frames = frame_scores.len/num_models;
   float mx_value;

   for (i=0; i<num_frames; i++) {

      // Find max value for current frame
      mx_value = frame_scores.data[i*num_models];
      i_mx = 0;
      for (j=1; j<num_models; j++) {
         if (frame_scores.data[i*num_models+j]> mx_value) {
            mx_value = frame_scores.data[i*num_models+j];
            i_mx = j;
         }
      }

      // initialize
      if (i==0) { 
         i_start = 0;
         curr_label = i_mx;
         continue;
      }

      // Check to see if label has changed
      if (i_mx != curr_label) {
         segments.push_back(Time_Interval(i_start*win_inc, i*win_inc, curr_label, mx_value));
         curr_label = i_mx;
         i_start = i;
      }

   }

   // Save the final interval
   segments.push_back(Time_Interval(i_start*win_inc, (num_frames-1)*win_inc, curr_label, mx_value));

} // find_intervals

static void smooth_frame_scores (vec<float> &frame_scores, int num_models, int win_len)
{
   vec<double> curr_output(num_models);
   vec<float> new_frame_scores(frame_scores.len);  // could do this in place, but easier logic this way

   int i, j, num_frames, win2;

   num_frames = frame_scores.len/num_models;
   win2 = win_len/2;

   if (num_frames < win_len) // don't smooth if there's not enough frames
      return;

   // Edge effect calculation -- use first win_len samples for averages
   for (i=0; i<num_models; i++)
      curr_output.data[i] = 0;
   for (i=0; i<win_len; i++) { 
      for (j=0; j<num_models; j++)
         curr_output.data[j] += frame_scores.data[i*num_models+j];
   }
   double mean = curr_output.mean();

   // Edge effect output
   for (i=0; i<win2; i++) {
      for (j=0; j<num_models; j++)
         new_frame_scores.data[i*num_models+j] = (float) ((curr_output.data[j]-mean)/((double) win_len));
   }

   // Standard output
   for (i=win2; i<(num_frames-win2); i++) {
      for (j=0; j<num_models; j++)
         curr_output.data[j] += frame_scores.data[(i+win2)*num_models+j]-frame_scores.data[(i-win2)*num_models+j];
      mean = curr_output.mean();
      for (j=0; j<num_models; j++)
         new_frame_scores.data[i*num_models+j] = (float) ((curr_output.data[j]-mean)/((double) win_len));
   }

   // Edge effects at the end
   for (i=num_frames-win2; i<num_frames; i++) {
      for (j=0; j<num_models; j++)
         new_frame_scores.data[i*num_models+j] = (float) ((curr_output.data[j]-mean)/((double) win_len));
   }

   // Copy result back into frame scores
   frame_scores = new_frame_scores;

} // smooth_frame_scores


