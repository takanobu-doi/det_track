#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#include <stdlib.h>
#include <math.h>
#include "para.h"

int hough_tra(unsigned int*** tpc_data, unsigned int*** hough_cnt);

void clear_hough_cnt(unsigned int*** hough_cnt);

int hough_tra(unsigned int*** tpc_data, unsigned int*** hough_cnt);

void find_hough_max(unsigned int*** hough_cnt,
		    unsigned int hough_max_pos[N_AC][2],
		    unsigned int hough_max_cnt[N_AC]);

int calc_pulse_integ_fadc(int*** fadc_data,
			  int pulse_integ[N_AC][N_FADC][FADC_MAX_PULSE],
			  int pulse_lead[N_AC][N_FADC][FADC_MAX_PULSE],
			  int pulse_width[N_AC][N_FADC][FADC_MAX_PULSE],
			  int integ_th);
