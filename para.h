#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#ifndef PARA
#define PARA

//#pragma once

const int N_TRY = 1;   // number of event generation for each (theta3,e3)
const int N_EVENT = 10000; //number of event generation
const int MAX_FILE_EVENT = 1000;

const int press = 500;   // for  500 hPa
//const int press = 1000;  // for 1000 hPa

const double FADC_WID    = 12.8;          // width of the 32 strpis (mm)
const double STRP_WID    =  0.4;          // width of the 1 strpi (mm)

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
const double STOP_X_MAX1 = FADC_WID*2.0;   // TPC boundary left (mm)
const double STOP_X_MIN2 = FADC_WID*5.0;   // TPC boundary right (mm)
const double STOP_X_MAX2 = STRP_WID*245.0; // TPC boundary right (mm)

const double STOP_Y_MIN  =   0.0;  // TPC boundary down (mm)
const double STOP_Y_MAX  = 110.0;  // TPC boundary up (mm)

const double STOP_Z_MIN  = 10*STRP_WID;       // TPC boundary upstream (mm)
const double STOP_Z_MAX  = 102.4-10*STRP_WID; // TPC boundary downsteram (mm)

const double beam_ene=750;
const double ene_10C = 780;
const double ene_12C = 1152;
const double AMU = 931.494095;
const double MassEx_10C = 15.6986;
const double MassEx_12C = 0;
const double MassEx_4He = 2.4249;
const double Mass_10C = AMU*10 + MassEx_10C;
const double Mass_12C = AMU*12 + MassEx_12C;
const double Mass_4He = AMU*4  + MassEx_4He;

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

const double GAS_GAIN = 1000.0;
const double TPC_THRESHOLD = 1.0;//default is 1.0

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
const double v_plate  = -4300.0;
const double y_grid   = 0.0;
const double v_grid   = -820.0;
const double E_FIELD  = (v_grid-v_plate)/(y_plate-y_grid);

const double step_size = 0.01; // step of avalanche

/* SRIM parameters */
const double Ratio_He   = 0.96;
const double Ratio_CO2  = 0.04;
const double Mass_He    =  4.0;
const double Mass_CO2   = 44.0;
const int    Charge_He  =  2;
const int    Charge_CO2 = 22;
const double density    = (1.1647e-4)*press/500.0;
const double W_Val      = 50.0; 
const double Fano_Factor= 1.0;
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

