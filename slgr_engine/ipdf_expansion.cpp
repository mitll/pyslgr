//
// Copyright 2012-2015 Massachusetts Institute of Technology. All rights
// reserved. Title to and copyright in the speaker, language, gender recognition
// software, any derivatives and any associated documentation is retained by
// MIT, subject to a non-exclusive royalty-free license to the U.S. Government
// to use the software for Government purposes as defined in Federal 
// Acquisition Regulation. Distribution to or use by third parties without 
// the prior written authorization of MIT is expressly prohibited.
//
// MIT MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED, WITH RESPECT
// TO THE FUNCTIONALITY OR USE OF THE SOFTWARE. ALL IMPLIED WARRANTIES, 
// INCLUDING BUT NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY OR FITNESS FOR 
// A PARTICULAR PURPOSE ARE EXPRESSLY EXCLUDED. FURTHER, MIT MAKES NO WARRANTY 
// THAT THE USE OF THE SOFTWARE COMPONENTS OR DOCUMENTATION WILL NOT INFRINGE 
// ANY PATENTS, COPYRIGHTS, TRADEMARKS, TRADE SECRETS OR OTHER RIGHTS OF ANY 
// PARTY.
//
// MIT DISCLAIMS ANY AND ALL LIABILITY FOR SPECIAL, INCIDENTAL, OR 
// CONSEQUENTIAL DAMAGES (INCLUDING LOSS OF PROFIT) ARISING OUT OF THE 
// INSTALLATION, USE, OPERATION OR SUPPORT OF THE SOFTWARE OR ASSOCIATED 
// DOCUMENTATION, EVEN IF MIT HAS BEEN APPRISED OF THE POSSIBILITY OF SUCH 
// DAMAGES.

//
// Standard header
//
#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

// Bill Campbell, 8/2008
// BC, Modified for IPDF, 12/2009
// Speech tools, 10/28/2010

#include "speech_tools.h"

//
// Classes
//
int IPDF_expansion::dim (void) 
{
	return exp_dim;
}

int IPDF_expansion::num_fea (void) 
{
	return exp_gmm.num_fea;
}

int IPDF_expansion::num_mix (void) 
{
	return exp_gmm.num_mix;
}

void IPDF_expansion::load_gmm_model (string model_file_name) {
	exp_gmm.load(model_file_name);
	exp_dim = exp_gmm.num_mix*exp_gmm.num_fea+1;
}

IPDF_expansion::IPDF_expansion (void) {
   exp_dim = 0;
	max_corank = 0;
}

IPDF_expansion::~IPDF_expansion (void) {
}

// Next two methods so IPDF_expansions can be used with STL
IPDF_expansion& IPDF_expansion::operator= (const IPDF_expansion &src) { 
   if (src.exp_dim!=0)
      throw ST_exception("Copying only allowed for uninitialized IPDF expansion.");
   exp_dim = 0;
   max_corank = 0;
   return *this;
}

IPDF_expansion::IPDF_expansion (const IPDF_expansion &src) {
   if (src.exp_dim!=0)
      throw ST_exception("Copying only allowed for uninitialized IPDF expansion.");
   exp_dim = 0;
   max_corank = 0;
}

vec<REAL_EXP> IPDF_expansion::expansion (Features &f, float rf) 
{
	if (!exp_gmm.is_loaded())
		throw ST_exception("Need to load GMM model.");

	vec<REAL_EXP> sum(f.num_outfeat()*exp_gmm.num_mix);
	vec<REAL_EXP> ec(exp_gmm.num_mix);
	vec<REAL_EXP> x_exp(exp_dim);

	// Compute suff stats
	exp_gmm.suff_stats_ns(f, sum, ec);

	// Unscaled expansion vector
	int i, j, k;
	x_exp.data[0] = 0;
	exp_gmm.map_adapt_mean(&(sum.data[0]), &(ec.data[0]), rf, &(x_exp.data[1])); 

	// Scaled expansion vector
	double tot_count = 0;
	double wgt;
	for (i=0; i<exp_gmm.num_mix; i++) {
		tot_count += (double) ec.data[i];
	}
	if (tot_count < 1.0e-8)
		tot_count = 1.0e-8;
	for (i=0; i<exp_gmm.num_mix; i++) {
		wgt = sqrt(ec.data[i]/tot_count);
		k = i*exp_gmm.num_fea;
		for (j=0; j<exp_gmm.num_fea; j++) {
			x_exp.data[k+j+1] -= exp_gmm.mean[k+j];
			x_exp.data[k+j+1] *= wgt*sqrt(exp_gmm.inv_cov[k+j]);
		}
	}

	return x_exp;

}

vec<REAL_EXP> IPDF_expansion::expansion_with_nap (Features &f, float rf, string nap_key, int corank) 
{
	if (!exp_gmm.is_loaded())
		throw ST_exception("Need to load GMM model.");

	vec<REAL_EXP> sum(f.num_outfeat()*exp_gmm.num_mix);
	vec<REAL_EXP> ec(exp_gmm.num_mix);
	vec<REAL_EXP> x_exp(exp_dim);

	// Compute suff stats
	exp_gmm.suff_stats_ns(f, sum, ec);

	// Unscaled expansion vector
	int i, j, k;
	x_exp.data[0] = 0;
	exp_gmm.map_adapt_mean(&(sum.data[0]), &(ec.data[0]), rf, &(x_exp.data[1])); 

	// Shift and scale
	double wgt;
	for (i=0; i<exp_gmm.num_mix; i++) {
		wgt = exp(0.5*exp_gmm.log_weight[i]);
		k = i*exp_gmm.num_fea;
		for (j=0; j<exp_gmm.num_fea; j++) {
			x_exp.data[k+j+1] -= exp_gmm.mean[k+j];
			x_exp.data[k+j+1] *= wgt*sqrt(exp_gmm.inv_cov[k+j]);
		}
	}

	// Compute NAP if necessary
	nap(x_exp, nap_key, corank);

	// Rescale
	scale_fixed_to_vm(x_exp, ec);

	return x_exp;

}

void IPDF_expansion::load_nap_projection (string projection_file_name, string key) {
	FILE *infile;

	if (!exp_gmm.is_loaded())
		throw ST_exception("Must load GMM model before projection is loaded");

	infile = fopen(projection_file_name.c_str(),"rb");
	if (infile==NULL) 
		throw ST_exception(string("unable to load NAP projection ")+projection_file_name);

	// Sanity check on size
	int st_size = file_len(projection_file_name.c_str());
	int num_entries = st_size/sizeof(REAL_EXP);
	int num_vec = num_entries/exp_dim;
	if (num_vec == 0) 
		throw ST_exception (string("empty NAP file ")+projection_file_name);
	if (sizeof(REAL_EXP)*num_vec*exp_dim != st_size) 
		throw ST_exception ("input file length is not a multiple of feature space dimension.");
	if ((max_corank!=0) && (max_corank!=num_vec))
		throw ST_exception("All NAP projections must have same max corank.");
	max_corank = num_vec;

	// Allocate space 
	vec<REAL_EXP> nap_proj(max_corank*exp_dim);

	// Now load in: headerless, raw doubles
	if (fread(nap_proj.data, sizeof(REAL_EXP), max_corank*exp_dim, infile)!=(max_corank*exp_dim))
		throw ST_exception("Unexpected EOF in NAP projection file.");
	fclose(infile);

	// Update class variables
	IPDF_expansion::max_corank = max_corank;
	nap_proj_map.insert(veckey::value_type(key, nap_proj));

}

void IPDF_expansion::nap (vec<REAL_EXP> &x_exp, string key, int corank) 
{
	int i, j;
	double sum;

	if (corank == 0)
		return;

	if (corank > max_corank)
		throw ST_exception("Corank is larger than max corank.");

	veckey::iterator it = nap_proj_map.find(key);
	if (it==nap_proj_map.end())
		throw ST_exception(string("NAP projection not found for the input key ")+key);
	vec<REAL_EXP> &nap_proj = it->second;

	for (i=0; i<corank; i++) {
		sum = 0;
		for (j=0; j<exp_dim; j++)
			sum += x_exp.data[j]*nap_proj.data[i*exp_dim+j];
		for (j=0; j<exp_dim; j++)
			x_exp.data[j] -= sum*nap_proj.data[i*exp_dim+j];
	}
}

void IPDF_expansion::scale_fixed_to_vm (vec<REAL_EXP> &x_exp, vec<REAL_EXP> &ec) {
	double wgt, tot_count;
	int i, j, k;

	tot_count = 0;
	for (i=0; i<exp_gmm.num_mix; i++) {
		tot_count += ec.data[i];
	}
	if (tot_count < 1.0e-8)
		tot_count = 1.0e-8;

	for (i=0; i<exp_gmm.num_mix; i++) {
		wgt = sqrt(ec.data[i]/tot_count)*exp(-0.5*exp_gmm.log_weight[i]);
		k = i*exp_gmm.num_fea;
		for (j=0; j<exp_gmm.num_fea; j++) {
			x_exp.data[k+j+1] *= wgt;
		}
	}
}

