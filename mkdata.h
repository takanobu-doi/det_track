#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#include <stdlib.h>
#include <math.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TSpline.h>

void clear_raw_wave(double*** raw_wave);

int get_upic_strp_num(double ele_end_pos);

int add_raw_wave(double ele_end_pos[4],   // input data 
		 double drift_time,       // input data
		 int ne,                  // input data
		 TGraph *wave,             // input wave template
		 double gain,             // input data
		 double*** raw_wave);     // output data

int add_raw_wave2(double ele_end_pos[4],   // input data 
		  double drift_time,       // input data
		  int ne,                  // input data
		  TSpline5 *wave_spline,   // input wave template
		  double gain,             // input data
		  double*** raw_wave);     // output data

void clear_fadc_data(int*** fadc_data);

int make_fadc_data(double*** raw_wave, int*** fadc_data);

void clear_tpc_data(unsigned int*** tpc_data);

int make_tpc_data(double*** raw_wave, unsigned int*** tpc_data, double threshold);

// original random generator for drift_electron()
double Uniform(void);
double rand_normal( double mu, double sigma );

// drift electron only take into account the diffusion
// use original random generator
int drift_electron(double start[4],
		   double driftv, double diff_tra, double diff_long,
		   double stop[4]);

// Use random generator from TRandom3 
int drift_electron2(double start[4],
		    double driftv, double diff_tra, double diff_long,
		    TRandom3 *gen_tra, TRandom3 *gen_long,
		    double stop[4]);

int add_fadc_noise(int*** fadc_data, TRandom3 *gen_noise, double noise);

double GetTrigTime(double vtx_y, double stop_y[], double driftv);
