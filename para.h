#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#ifndef _PARA_H_
#define _PARA_H_

#include <map>
#include <string>
//#pragma once

const int N_TRY = 1;   // number of event generation for each (theta3,e3)
const int N_EVENT = 1000; //number of event generation
const int MAX_FILE_EVENT = 1000;

const char detection_gas[128] = "CH4-1"; // CH4 iso-C4H10 CH4-* iso-C4H10-*
const int press = 500;   // hPa
//const int press = 1000;  // for 1000 hPa
const std::string FILENAME = "test";

const Int_t nAlpha = 2; //number of decay products

const double FADC_WID    = 12.8;          // width of the 32 strips (mm)
const double STRP_WID    =  0.4;          // width of the 1 strip (mm)

const double THETA3_START = 0.0;   // minimum LAB angle (degree)
const double THETA3_STOP  = 100.1; // maximum LAB angle (degree)
const double THETA3_STEP  = 1.0;   // LAB angle step (degree)
const int N_THETA3 = 1000;
const int THETA3_BIN = (THETA3_STOP-THETA3_START)/THETA3_STEP;

const double E3_START = 0.0;       // minimum LAB energy (MeV)
const double E3_STOP  = 3.0;      // maximum LAB energy (MeV)
const double E3_STEP  = 0.1;      // LAB energy step (MeV)
const int N_E3 = 100;
const int E3_BIN = (E3_STOP-E3_START)/E3_STEP;

const double PHI3_START = 0.0;     // minimum phi angle (degree)
const double PHI3_STOP  = 360.0;   // maximum phi angle (degree)
const double PHI3_STEP  = 5.0;    // maximum phi angle (degree)

const double RANGE_RESO = 2.0;     // range resolution in sigma (mm)

const double VTX_X_MEAN  = 102.4/2.0;  // beam center in horizontal direction (mm)
const double VTX_X_SIGMA =   0.0;      // beam spread in horizontal direction (mm)

const double VTX_Y_MEAN  =  55.0;      // beam center in vertical direction (mm)
const double VTX_Y_SIGMA =   0.0;      // beam spread in vertical direction (mm)

const double VTX_Z_START = FADC_WID*1.0;  // minimum vertex position (mm)
const double VTX_Z_STOP  = FADC_WID*7.0;  // maximum vertex position (mm)
const double VTX_Z_STEP  = 0.1;           // vertex position step (mm)

const double STOP_X_MIN1 = STRP_WID*10.0;  // TPC boundary left (mm)
//****const double STOP_X_MAX1 = FADC_WID*2.0;   // TPC boundary left (mm) <-- FADC_WID*3.0 ?
//****const double STOP_X_MIN2 = FADC_WID*5.0;   // TPC boundary right (mm)
const double STOP_X_MAX1 = FADC_WID*4.0;  // changed by T.Doi 2019/3/31
const double STOP_X_MIN2 = FADC_WID*4.0;  // changed by T.Doi 2019/3/31
const double STOP_X_MAX2 = STRP_WID*245.0; // TPC boundary right (mm)

const double STOP_Y_MIN  =   0.0;  // TPC boundary down (mm)
const double STOP_Y_MAX  = 110.0;  // TPC boundary up (mm)

const double STOP_Z_MIN  = 10*STRP_WID;       // TPC boundary upstream (mm)  take 10 strips mergin
const double STOP_Z_MAX  = 102.4-10*STRP_WID; // TPC boundary downsteram (mm)

const double beam_ene = 750.;
const double ene_10C = 780;
const double ene_12C = 1152;
const double AMU = 931.494095;
const double MassEx_10C = 15.6986;
const double MassEx_12C = 0;
const double MassEx_4He = 2.4249;
const double MassEx_n = 8.0713;
const double MassEx_p = 7.2889;
const double Mass_10C = AMU*10 + MassEx_10C;
const double Mass_12C = AMU*12 + MassEx_12C;
const double Mass_4He = AMU*4  + MassEx_4He;
const double Mass_n = AMU*1 + MassEx_n;
const double Mass_p = AMU*1 + MassEx_p;

const int max_particle = 3; // 0:4He 1:10C 2:12C 3:16O

const int N_AC = 2;
const int ANODE = 0;
const int CATHODE = 1;
const int N_TCLK = 1024;
const int N_STRP = 256;
const int N_FADC = 8;
const int FADC_MAX_PULSE = 4;
const int N_ACLK =256;
const double TPC_SIZE = 102.4;
const double mVToFADC = 256/800.0;
const double FADC_NOISE = 1.5;
const int FADC_PULSE_TH = 2;

//****const double GAS_GAIN = 1000.0; // this is default value
const double GAS_GAIN = 2000.0;
const double TPC_THRESHOLD = 1.;//default is 1.0
const double SAMPLING_RATIO = 0.1; // GHz default is 0.1 GHz
const double DELAY_RATIO = 90.; // % delay time of read out vs time window
const double BUFF_TIME = N_TCLK/SAMPLING_RATIO*(1.-DELAY_RATIO/100.); // ns buffer time of trig

///////////////////////////
//  garfield parameters  //
///////////////////////////
const double mmTocm = 0.1;
const double cmTomm = 10.0;
const double center_x = 10.24/2.0;
const double center_y = 7.0;
const double center_z = 10.24/2.0;
const double half_x   = 15.0;
const double half_y   = 15.0;
const double half_z   = 15.0;
const double y_plate  = 14.0;
//const double v_plate  = -4300.0; this is default value in furuno exp.
const double v_plate = -6300.0;
const double y_grid   = 0.0;
const double v_grid   = -820.0;
const double E_FIELD  = (v_grid-v_plate)/(y_plate-y_grid);

const double step_size = 0.001; // step of avalanche

#ifdef SRIM_PARA
/* SRIM parameters */
//****const double Ratio_He   = 0.96;
//****const double Ratio_CO2  = 0.04;
std::map<const std::string,const double> Mass = {
//****  {"He",4.0},
//****  {"CO2",44.0},
  {"iso-C4H10",58.}, // added by T. Doi 2019/3/28
  {"CH4",16.},
  {"iso-C4H10-90",58.*.90+1.*(1.-.90)},
  {"CH4-90",16*.90+1.*(1.-.90)},
  {"iso-C4H10-80",58.*.80+1.*(1.-.80)},
  {"CH4-80",16*.80+1.*(1.-.80)},
  {"iso-C4H10-70",58.*.70+1.*(1.-.70)},
  {"CH4-70",16*.70+1.*(1.-.70)},
  {"iso-C4H10-60",58.*.60+1.*(1.-.60)},
  {"CH4-60",16*.60+1.*(1.-.60)},
  {"iso-C4H10-50",58.*.50+1.*(1.-.50)},
  {"CH4-50",16*.50+1.*(1.-.50)},
  {"iso-C4H10-40",58.*.40+1.*(1.-.40)},
  {"CH4-40",16*.40+1.*(1.-.40)},
  {"iso-C4H10-30",58.*.30+1.*(1.-.30)},
  {"CH4-30",16*.30+1.*(1.-.30)},
  {"iso-C4H10-20",58.*.20+1.*(1.-.20)},
  {"CH4-20",16*.20+1.*(1.-.20)},
  {"iso-C4H10-10",58.*.10+1.*(1.-.10)},
  {"CH4-10",16*.10+1.*(1.-.10)},
  {"iso-C4H10-5",58.*.05+1.*(1.-.05)},
  {"CH4-5",16.*.05+1.*(1.-.05)},
  {"iso-C4H10-1",58.*.01+1.*(1.-.01)},
  {"CH4-1",16.*.01+1.*(1.-.01)}
};

std::map<const std::string,const double> Charge = {
//****  {"He",2},
//****  {"CO2",22},
  {"iso-C4H10",34}, //added by T. Doi 2019/3/28
  {"CH4",10},
  {"iso-C4H10-90",34*.90+1.*(1.-.90)},
  {"CH4-90",10*.90+1.*(1.-.90)},
  {"iso-C4H10-80",34*.80+1.*(1.-.80)},
  {"CH4-80",10*.80+1.*(1.-.80)},
  {"iso-C4H10-70",34*.70+1.*(1.-.70)},
  {"CH4-70",10*.70+1.*(1.-.70)},
  {"iso-C4H10-60",34*.60+1.*(1.-.60)},
  {"CH4-60",10*.60+1.*(1.-.60)},
  {"iso-C4H10-50",34*.50+1.*(1.-.50)},
  {"CH4-50",10*.50+1.*(1.-.50)},
  {"iso-C4H10-40",34*.40+1.*(1.-.40)},
  {"CH4-40",10*.40+1.*(1.-.40)},
  {"iso-C4H10-30",34*.30+1.*(1.-.30)},
  {"CH4-30",10*.30+1.*(1.-.30)},
  {"iso-C4H10-20",34*.20+1.*(1.-.20)},
  {"CH4-20",10*.20+1.*(1.-.20)},
  {"iso-C4H10-10",34*.10+1.*(1.-.10)},
  {"CH4-10",10*.10+1.*(1.-.10)},
  {"iso-C4H10-5",34*.05+1.*(1.-.05)},
  {"CH4-5",10*.05+1.*(1.-.05)},
  {"iso-C4H10-1",34*.01+1.*(1.-.01)},
  {"CH4-1",10*.01+1.*(1.-.01)}
};
//****const double density    = (1.1647e-4)*press/500.;
std::map<const std::string,const double> density = {
  {"iso-C4H10",0.0023846*press/1000.}, // g/cm3 changed by T. Doi 2019/3/28
  {"CH4",0.00065819*press/1000.},
  {"iso-C4H10-90",(0.0023846*.90+0.00008271*(1.-.90))*press/1000.},
  {"CH4-90",(0.00065819*.90+0.00008271*(1.-.90))*press/1000.},
  {"iso-C4H10-80",(0.0023846*.80+0.00008271*(1.-.80))*press/1000.},
  {"CH4-80",(0.00065819*.80+0.00008271*(1.-.80))*press/1000.},
  {"iso-C4H10-70",(0.0023846*.70+0.00008271*(1.-.70))*press/1000.},
  {"CH4-70",(0.00065819*.70+0.00008271*(1.-.70))*press/1000.},
  {"iso-C4H10-60",(0.0023846*.60+0.00008271*(1.-.60))*press/1000.},
  {"CH4-60",(0.00065819*.60+0.00008271*(1.-.60))*press/1000.},
  {"iso-C4H10-50",(0.0023846*.50+0.00008271*(1.-.50))*press/1000.},
  {"CH4-50",(0.00065819*.50+0.00008271*(1.-.50))*press/1000.},
  {"iso-C4H10-40",(0.0023846*.40+0.00008271*(1.-.40))*press/1000.},
  {"CH4-40",(0.00065819*.40+0.00008271*(1.-.40))*press/1000.},
  {"iso-C4H10-30",(0.0023846*.30+0.00008271*(1.-.30))*press/1000.},
  {"CH4-30",(0.00065819*.30+0.00008271*(1.-.30))*press/1000.},
  {"iso-C4H10-20",(0.0023846*.20+0.00008271*(1.-.20))*press/1000.},
  {"CH4-20",(0.00065819*.20+0.00008271*(1.-.20))*press/1000.},
  {"iso-C4H10-10",(0.0023846*.10+0.00008271*(1.-.10))*press/1000.},
  {"CH4-10",(0.00065819*.10+0.00008271*(1.-.10))*press/1000.},
  {"iso-C4H10-5",(0.0023846*.05+0.00008271*(1.-.05))*press/1000.},
  {"CH4-5",(0.00065819*.05+0.00008271*(1.-.05))*press/1000.},
  {"iso-C4H10-1",(0.0023846*.01+0.00008271*(1.-.01))*press/1000.},
  {"CH4-1",(0.00065819*.01+0.00008271*(1.-.01))*press/1000.}
};
//const double W_Val      = 50.0;
std::map<const std::string,const double> W_Val = {
  {"iso-C4H10",10.68},
  {"CH4",12.61},
  {"iso-C4H10-90",10.68*.90+13.598*(1.-.90)},
  {"CH4-90",12.61*.90+13.598*(1.-.90)},
  {"iso-C4H10-80",10.68*.80+13.598*(1.-.80)},
  {"CH4-80",12.61*.80+13.598*(1.-.80)},
  {"iso-C4H10-70",10.68*.70+13.598*(1.-.70)},
  {"CH4-70",12.61*.70+13.598*(1.-.70)},
  {"iso-C4H10-60",10.68*.60+13.598*(1.-.60)},
  {"CH4-60",12.61*.60+13.598*(1.-.60)},
  {"iso-C4H10-50",10.68*.50+13.598*(1.-.50)},
  {"CH4-50",12.61*.50+13.598*(1.-.50)},
  {"iso-C4H10-40",10.68*.40+13.598*(1.-.40)},
  {"CH4-40",12.61*.40+13.598*(1.-.40)},
  {"iso-C4H10-30",10.68*.30+13.598*(1.-.30)},
  {"CH4-30",12.61*.30+13.598*(1.-.30)},
  {"iso-C4H10-20",10.68*.20+13.598*(1.-.20)},
  {"CH4-20",12.61*.20+13.598*(1.-.20)},
  {"iso-C4H10-10",10.68*.10+13.598*(1.-.10)},
//****  {"iso-C4H10-10",10.},
  {"CH4-10",12.61*.10+13.598*(1.-.10)},
  {"iso-C4H10-5",10.68*.05+13.598*(1.-.05)},
  {"CH4-5",12.61*.05+13.598*(1.-.05)},
  {"iso-C4H10-1",10.68*.01+13.598*(1.-.01)},
  {"CH4-1",12.61*.01+13.598*(1.-.01)}
};
#endif

const double Fano_Factor= 1.0;
//****const int Cluster_Size = 100;
const int Cluster_Size = 100;
const double eVToMeV = 1.0e-6;
const double MeVToeV = 1.0e6;

/* analysis parameter */
const int DIV_HOUGH_X = 180;
const int DIV_HOUGH_Y = 256;
const double PI = 3.1415926535;
const int HOUGH_TH_A = 160;
const int HOUGH_TH_C = 120;
const double HOUGH_RECOIL_ANG_A = 12.0;

/* original ionization */
const double ion_step = 0.1;   // in mm
const double eneloss_fac_alpha = 1.1647e-2;


#endif
