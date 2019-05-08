/* standard header files */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <dirent.h>
#include <cstring>
#include <string.h>
#include <unistd.h>

/* header files for ROOT */
#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TMath.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSpline.h>

/* header files for garfield++ */
#include "MediumMagboltz.hh"
#include "FundamentalConstants.hh"
#include "ComponentAnalyticField.hh"
#include "GeometrySimple.hh"
#include "SolidBox.hh"
#include "Sensor.hh"
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "ViewField.hh"
#include "ViewCell.hh"
#include "ViewDrift.hh"
#include "DriftLineRKF.hh"
#include "TrackSrim.hh"
using namespace Garfield;

#define VIEW_EVE 1
#define ALPHA_TRACK 1
#define BEAM_TRACK 1
#define ORIGINAL_DRIFT 1
//#define ONE_ELE 1

#include "kinema.h"
#include "mkdata.h"
#include "anatpc.h"
#include "para.h"
#include "database.hpp"
#include "kinema_doi.hpp"
#include "nuclear.hpp"

int judge_stop_inside(double stop_x, double stop_y, double stop_z);

int main(int argc, char *argv[]){
  int i,j,k;

  /* get the working directory */
//**  ******** changed by T.Doi 2019/3/4 **********
//**  const char* det_track = "/det_track";
//**  static char exename[1024];
//**  char workdir[1024];
//**  readlink( "/proc/self/exe", exename, sizeof(exename)-1 );
//**  printf("exe=%s\n", exename);
//**  strncpy(workdir, exename, strlen(exename)-strlen(det_track));
//**  printf("workdir=%s\n", workdir);
//**
  char workdir[1024];
  getcwd(workdir,1024);
  
  int index=1;
  char outfname[512];

  int theta_start, theta_stop;
  if(argc!=2){
    index=1;
    theta_start = THETA3_START;
    theta_stop = THETA3_STOP;    
  }
  if(argc==2){
    index = atoi(argv[1]);
    theta_start = ((int)(index*THETA3_STEP))%100;
    theta_stop = theta_start + THETA3_STEP;
  }

  printf("index=%d\n", index);
  sprintf(outfname, "%s/root/out%d.root", workdir, index);
  if(argc==1) sprintf(outfname, "%s/out.root", workdir);
  printf("outfile: %s\n", outfname);
  
  int total_loop=1;
  //  total_loop*=(int)((theta_stop-theta_start)/THETA3_STEP);
  total_loop*=(int)((E3_STOP-E3_START)/E3_STEP);
  total_loop*=(int)((PHI3_STOP-PHI3_START)/PHI3_STEP);
  total_loop*=N_TRY;
  
  printf("total loop= %d\n", total_loop);
    
  /* initialize random seed (always same init) */
  TRandom3 *rndm = new TRandom3();
  rndm->SetSeed(index);
  
  /* read energy to range table */
  char grfilename[512], gr2filename[512];
  sprintf(grfilename,  "%s/tables/ene_to_range_%d.dat", workdir, press);
  sprintf(gr2filename, "%s/tables/range_to_ene_%d.dat", workdir, press);  
  printf("range table: %s\n", grfilename);
  TGraph *gr_range;
  gr_range = new TGraph(grfilename);
  TGraph *gr_ene;
  gr_ene = new TGraph(gr2filename);
  
  /* read energy loss table */
  char enelossfname_alpha[512];
  sprintf(enelossfname_alpha, "%s/tables/eneloss_alpha_%d.dat", workdir, press);
  TGraph *gr_eneloss_alpha = new TGraph(enelossfname_alpha);
  

  double theta3_deg;   // recoil angle in LAB (degree)
  double theta3_rad;   // recoil angle in LAB (radian)  
  double e3;           // recoil energy in LAB (MeV)
  double e3_rec;       // reconstructed recoil energy with finite resolution in LAB (MeV)  
  double range;        // recoil range (mm)
  double range_rec;    // reconstructed recoil range with finite resolution (mm)  
  double phi3_deg;     // recoil phi angle (degree)  
  double phi3_rad;     // recoil phi angle (radian)  
  double vtx_x,vtx_y,vtx_z;     // vertex position (mm)
  double dx,dy,dz;              // position difference between vertex and stop points
  double dr;
  double stop_x,stop_y,stop_z;  // recoil particle stop position (mm)
  double ex;                    // excitation energy of 10C
  double theta_cm;              // angle in center of mass

  double beam_start_x, beam_start_y, beam_start_z;
  
  int itheta3, ie3;
  int itry;
  int stop_inside;

  int total_cnt;    // total number of event for each (theta,e3)
  int stop_cnt;     // number of stop inside TPC event for each (theta,e3)
  double stop_eff;  // ratio of the stop_cnt/total_cnt
  int neve;         // number of total generated event

  kinema a;
  a.setparticles("10c","12c","12c","10c");

  /* garfield parameters */
  double cluster_pos[4];
  double ele_end_pos[4];  
  double drift_time;
  int drift_status;
  int n_cluster;
  int tot_ne;
  double t_0=0;
  int ne;
  double xc, yc, zc, tc, ec, ekin;  

  double recoil_track_a_x[2], recoil_track_a_y[2];
  double recoil_track_c_x[2], recoil_track_c_y[2];  
  
  /* Hough analysis parameters */
  int pulse_integ[N_AC][N_FADC][FADC_MAX_PULSE];
  int pulse_lead[N_AC][N_FADC][FADC_MAX_PULSE];
  int pulse_width[N_AC][N_FADC][FADC_MAX_PULSE];  
  unsigned int hough_max_pos[N_AC][2];  
  unsigned int hough_max_cnt[N_AC];
  double hough_theta[N_AC];
  double hough_r[N_AC];
  
  double line_para[N_AC][2];
  double line_y[N_AC][2];
  
  unsigned int hough_scatt_flag;
  
//**  TH2D *h100 = new TH2D("h100", "efficiency",
//**			THETA3_BIN, THETA3_START, THETA3_STOP,
//**			E3_BIN, E3_START, E3_STOP);
//**
//**  TH2D *h101 = new TH2D("h101", "all: e3 vs theta3",
//**			THETA3_BIN, THETA3_START, THETA3_STOP,
//**			E3_BIN, E3_START, E3_STOP);
//**
//**  TH2D *h102 = new TH2D("h102", "stop:e3 vs theta3",
//**			THETA3_BIN, THETA3_START, THETA3_STOP,
//**			E3_BIN, E3_START, E3_STOP);
//**
//**  TH2D *h103 = new TH2D("h103", "all: event theta_cm vs ex",
//**			100,-5,45, 20,0,10);
//**  TH2D *h104 = new TH2D("h104", "stop: event theta_cm vs ex",
//**			100,-5,45, 20,0,10);
//**  
//**  TH2D *h105 = new TH2D("h105", "all: event e3 vs range",
//**			300,0,150, 200,0,5);
//**  TH2D *h106 = new TH2D("h106", "stop: event e3 vs range",
//**			300,0,150, 200,0,5);
//**
//**  TH1D *h200 = new TH1D("h200", "vertex z gated", 100,0,100);
//**  
//**  TH1F *h201 = new TH1F("h201", "range, stopped", 200,0,200);
//**  TH1F *h202 = new TH1F("h202", "phi, stopped", 360,0,360);
//**  TH2F *h203 = new TH2F("h203", "vertical length vs phi, stopped",
//**			360,0,360, 200,-100,100);    
//**  TH1F *h204 = new TH1F("h204", "e3, stopped", 500,0,10);
//**
//**
//**  TH2F *h300 = new TH2F("h300", "stop && Hough:e3 vs theta3",
//**			THETA3_BIN, THETA3_START, THETA3_STOP,
//**			E3_BIN, E3_START, E3_STOP);
//**
//**  TH2F *h301 = new TH2F("h301", "stop && Hough:theta_cm vs ex",
//**			100,-5,45, 20,0,10);
//**
//**  TH1F *h302 = new TH1F("h302", "stop && Hough:phi", 360,0,360);
//**  TH1D *h303 = new TH1D("h303", "stop && Hough:vertex", 100,0,100);
//**  
//**  /* electron start points */
//**  TH2F *h1000 = new TH2F("h1000", "cluster position Anode",   1024,0,102.4, 1000,0,150);
//**  TH2F *h1001 = new TH2F("h1001", "cluster position Cathode", 1024,0,102.4, 1000,0,150);  
//**  TH1F *h1002 = new TH1F("h1002", "cluster position Z", 256,0,102.4);
//**  TH1F *h1003 = new TH1F("h1003", "cluster position X", 256,0,102.4);  
//**
//**  /* electron stop points */
//**  TH2F *h1010 = new TH2F("h1010", "drift end position Anode",  1024,0,102.4, 1000,0,150);
//**  TH2F *h1011 = new TH2F("h1011", "drift end position Cathode",1024,0,102.4, 1000,0,150);  
//**  TH1F *h1012 = new TH1F("h1012", "drift end position Z", 256,0,102.4);
//**  TH1F *h1013 = new TH1F("h1013", "drift end position X", 256,0,102.4);  
//**  
//**  TH2D *h1100 = new TH2D("h1100", "raw wave Anode", 256,0,256, 1024,0,1024);
//**  TH2D *h1101 = new TH2D("h1101", "raw wave Cathode", 256,0,256, 1024,0,1024);  
//**
//**  TH2D *h1200 = new TH2D("h1200", "FADC Anode", 8,0,8, 256,0,256);
//**  TH2D *h1201 = new TH2D("h1201", "FADC wave Cathode", 8,0,8, 256,0,256);  
//**
//**  TH1F hfadc[N_AC][N_FADC];
//**  for(int i=0; i<N_AC; i++){
//**    for(int j=0; j<N_FADC; j++){
//**      hfadc[i][j] = TH1F(Form("hfadc%d%d",i,j), Form("wave form ch [%d][%d]",i,j),
//**                         256,0,256);
//**    }
//**  }
//**  
//**  TH2D *h1300 = new TH2D("h1300", "TPC data Anode", 256,0,256, 1024,0,1024);
//**  TH2D *h1301 = new TH2D("h1301", "TPC data Cathode", 256,0,256, 1024,0,1024);  
//**
//**  TH2I *h2000 = new TH2I("h2000", "Hough space Anode",
//**			 DIV_HOUGH_X,0,360, DIV_HOUGH_Y,0,1024);
//**  TH2I *h2001 = new TH2I("h2001", "Hough space Cathode",
//**			 DIV_HOUGH_X,0,360, DIV_HOUGH_Y,0,1024);
//**  
//**#ifdef VIEW_EVE
//**  TApplication app("app", &argc, argv);
//**  TStyle *style = new TStyle("Plain", "style1");
//**  style->SetPalette(55);
//**  style->SetOptStat(0);
//**  style->cd();
//**
//**  TCanvas *c1 = new TCanvas("name1", "2-dim cluster distribution", 900,600);
//**  c1->SetLeftMargin(0.13);
//**  c1->SetFillStyle(0);
//**  TCanvas *c2 = new TCanvas("name2", "1-dim cluster distribution", 900,600);
//**  c2->SetLeftMargin(0.13);
//**  c2->SetFillStyle(0);
//**
//**  TCanvas *c3 = new TCanvas("name3", "2-dim electron end distribution", 900,600);
//**  c3->SetLeftMargin(0.13);
//**  c3->SetFillStyle(0);
//**  TCanvas *c4 = new TCanvas("name4", "1-dim electron end distribution", 900,600);
//**  c4->SetLeftMargin(0.13);
//**  c4->SetFillStyle(0);
//**
//**  TCanvas *c5 = new TCanvas("name5", "2-dim wave", 900,600);
//**  c5->SetLeftMargin(0.13);
//**  c5->SetFillStyle(0);
//**
//**  TCanvas *c6 = new TCanvas("name6", "TPC track", 900,600);
//**  c6->SetLeftMargin(0.13);
//**  c6->SetFillStyle(0);
//**
//**  TCanvas *c10 = new TCanvas("name10", "Hough space", 900,600);
//**  c10->SetLeftMargin(0.13);
//**  c10->SetFillStyle(0);
//**#endif
  
  ///////////////////////////
  //   garfield settings   //
  ///////////////////////////

  /* Read Magboltz file */
  char magfname[512];
  sprintf(magfname, "%s/tables/He-96_CO2-4_%d.gas", workdir, press);

  MediumMagboltz* gas = new MediumMagboltz();
  int magres = gas->LoadGasFile(magfname);
  if(magres==0){
    printf("error in reading Magboltz file: %s\n", magfname);
    exit(0);
  } 
  
  /* get gas property */
  std::vector<double> efields, bfields, angles;
  gas->GetFieldGrid(efields, bfields, angles);
  double ef[efields.size()];
  double velocity[efields.size()];
  double diff_t[efields.size()], diff_l[efields.size()];
  double town_a[efields.size()], attach[efields.size()];
  for(i=0; i<efields.size(); i++){
    gas->GetElectronVelocityE(i, 0, 0, velocity[i]);
    gas->GetElectronTransverseDiffusion(i, 0, 0, diff_t[i]);
    gas->GetElectronLongitudinalDiffusion(i, 0, 0, diff_l[i]);
    gas->GetElectronTownsend(i, 0, 0, town_a[i]);
    gas->GetElectronAttachment(i, 0, 0, attach[i]);
    
    ef[i]= efields[i];
    velocity[i]=velocity[i]*1000;
    diff_t[i]=diff_t[i]*10000;
    diff_l[i]=diff_l[i]*10000;
    
    if(town_a[i]<0) town_a[i]=0;
    if(attach[i]<0) attach[i]=0;
  }  

  /* calculate drift velocity */
  TGraph *gr_driftv = new TGraph(efields.size(), ef, velocity);
  double driftv = gr_driftv->Eval(E_FIELD)*0.1;
  
  /* calculate trans diffusion */ 
  TGraph *gr_diff_tra = new TGraph(efields.size(), ef, diff_t);
  double diff_tra = gr_diff_tra->Eval(E_FIELD);

  /* calculate long diffusion */ 
  TGraph *gr_diff_long = new TGraph(efields.size(), ef, diff_l);
  double diff_long = gr_diff_long->Eval(E_FIELD);

  printf("E: %.1f V/cm, driftv: %.3f, diff_tra: %.1f, long_diff: %.1f\n",
	 E_FIELD, driftv, diff_tra, diff_long);
  
  /* Define the gas volume */
  GeometrySimple* geo = new GeometrySimple();
  SolidBox* box = new SolidBox(center_x, center_y, center_z,
                               half_x, half_y, half_z);
  geo->AddSolid(box, gas);  

  /* Setup the electric field */
  ComponentAnalyticField* comp = new ComponentAnalyticField();

  /* TPC cage */
  comp->AddPlaneY(y_plate, v_plate, "plate");
  comp->AddPlaneY(y_grid,  v_grid,  "grid" );
  comp->SetGeometry(geo);

  /* Create a sensor */
  Sensor* sensor = new Sensor();
  sensor->AddComponent(comp);
  sensor->SetArea(center_x-half_x, center_y-half_y, center_z-half_z,
		  center_x+half_x, center_y+half_y, center_z+half_z);

  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  //  drift->SetDistanceSteps(step_size);
  
  /* Read Srim file */
  char srimfname_alpha[512], srimfname_10c[512];
  sprintf(srimfname_alpha, "%s/tables/12C_HeCO2_96_4_%d.srim", workdir, press);
  sprintf(srimfname_10c,   "%s/tables/10C_HeCO2_96_4_%d.srim", workdir, press);  

  TrackSrim* srim_alpha = new TrackSrim();
  srim_alpha->SetSensor(sensor);
  int srimres = srim_alpha->ReadFile(srimfname_alpha);
  if(srimres==0){
    printf("error in reading SRIM file: %s\n", srimfname_alpha);
    exit(0);
  }

  srim_alpha->SetWorkFunction(W_Val);  // in eV
  srim_alpha->SetFanoFactor(Fano_Factor);
  srim_alpha->SetModel(4);
  srim_alpha->SetAtomicMassNumbers(Ratio_He*Mass_He + Ratio_CO2*Mass_CO2,
			     Ratio_He*Charge_He + Ratio_CO2*Charge_CO2);
  srim_alpha->SetDensity(density);

  srim_alpha->SetTargetClusterSize(Cluster_Size);  
  srim_alpha->DisableTransverseStraggling();
  srim_alpha->DisableLongitudinalStraggling();

  printf("************************************************************\n");
  srim_alpha->Print();
  printf("************************************************************\n");

  /* read srim file for 10C beam */
#ifdef BEAM_TRACK
  TrackSrim* srim_10c = new TrackSrim();
  srim_10c->SetSensor(sensor);
  int srimres2 = srim_10c->ReadFile(srimfname_10c);
  if(srimres2 ==0){
    printf("error in reading SRIM file: %s\n", srimfname_10c);
    exit(0);
  }
  srim_10c->SetWorkFunction(W_Val);  // in eV
  srim_10c->SetFanoFactor(Fano_Factor);
  //  srim_10c->SetModel(3);
  srim_10c->SetAtomicMassNumbers(Ratio_He*Mass_He + Ratio_CO2*Mass_CO2,
				 Ratio_He*Charge_He + Ratio_CO2*Charge_CO2);
  srim_10c->SetDensity(density);

  srim_10c->SetTargetClusterSize(1);  
  srim_10c->DisableTransverseStraggling();
  srim_10c->DisableLongitudinalStraggling();

  printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  srim_10c->Print();
  printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  srim_10c->SetKineticEnergy((beam_ene+0.0001)*1.0e6);
#endif
  
  double beam_dx=0;
  double beam_dy=0;
  double beam_dz=-1.0;

  // random geenrators
  /* diffusion random */
  TRandom3 *gen_tra  = new TRandom3();
  TRandom3 *gen_long = new TRandom3();  
  gen_tra->SetSeed(time(NULL));
  gen_tra->SetSeed(time(NULL)+1);

  /* noise */
  TRandom3 *gen_fadc_noise = new TRandom3();
  gen_fadc_noise->SetSeed(time(NULL));

  srand((unsigned)time(NULL));
  
  //////////////////////////////////////////
  //   prepare for TPC data conversion    //
  //////////////////////////////////////////
  
  /* analog pulse data for 256 strips */
  double*** raw_wave;
  raw_wave = (double***)malloc(sizeof(double**)*N_AC);
  if(raw_wave==NULL){
    printf("error in malloc of raw_wave");
    exit(1);
  }
  for(i=0; i<N_AC; i++){
    raw_wave[i] = (double**)malloc(sizeof(double*)*N_TCLK);
    if(raw_wave[i]==NULL){
      printf("error in malloc of raw_wave");
      exit(1);
    }
    for(j=0; j<N_TCLK; j++){
      raw_wave[i][j] = (double*)malloc(sizeof(double)*N_STRP);
      if(raw_wave[i][j]==NULL){
	printf("error in malloc of raw_wave");
	exit(1);
      }
    }
  }
  clear_raw_wave(raw_wave);

  /* analog pulse data for FADC */
  int*** fadc_data;
  fadc_data = (int***)malloc(sizeof(int**)*N_AC);
  if(fadc_data==NULL){
    printf("error in malloc of fadc_data\n");
    exit(1);
  }
  for(i=0; i<N_AC; i++){
    fadc_data[i] = (int**)malloc(sizeof(int*)*N_FADC);
    if(fadc_data[i]==NULL){
      printf("error in malloc of fadc_data\n");
      exit(1);
    }
    for(j=0; j<N_FADC; j++){
      fadc_data[i][j] = (int*)malloc(sizeof(int)*N_ACLK);
      if(fadc_data[i][j]==NULL){
	printf("error in malloc of fadc_data\n");
	exit(1);
      }
    }
  }
  clear_fadc_data(fadc_data);
  
  /* TPC hit data */
  unsigned int*** tpc_data;
  tpc_data = (unsigned int***)malloc(sizeof(unsigned int**)*N_AC);
  if(tpc_data==NULL){
    printf("error in malloc of tpc_data\n");
    exit(1);
  }
  for(i=0; i<N_AC; i++){
    tpc_data[i] = (unsigned int**)malloc(sizeof(unsigned int*)*N_TCLK);
    if(tpc_data[i]==NULL){
      printf("error in malloc of tpc_data\n");
      exit(1);
    }
    for(j=0; j<N_TCLK; j++){
      tpc_data[i][j] = (unsigned int*)malloc(sizeof(unsigned int)*N_STRP);
      if(tpc_data[i][j]==NULL){
	printf("error in malloc of tpc_data\n");
	exit(1);
      }
    }      
  }
  clear_tpc_data(tpc_data);  

  /* wave template data */
  char wavefname[512];
  sprintf(wavefname, "%s/tables/wave_temp.dat", workdir);
  FILE* wave;
  wave = fopen(wavefname, "r");
  if(wave == NULL){
    printf("error in reading wave file: %s\n", wavefname);
    exit(0);
  }
  fclose(wave);
  TGraph *wave_temp = new TGraph(wavefname);
  wave_temp->SetBit(1);
  
  TSpline5* wave_spline=new TSpline5("wave_spline", wave_temp);

  const int filesize = N_EVENT/MAX_FILE_EVENT;
  FILE *fp_t[filesize];
  FILE *fp_r[filesize];
  char filename[1024];
  for(int ii=0;ii<filesize;ii++){
    sprintf(filename,"12C_track-%d.dat",ii);
    fp_t[ii] = fopen(filename,"w");
    sprintf(filename,"12C_result-%d.dat",ii);
    fp_r[ii] = fopen(filename,"w");
  }

  //////////////////////////////////////////
  //       prepare for TPC analysis       //
  //////////////////////////////////////////

//**  /* Hough space data */
//**  unsigned int*** hough_cnt;
//**  hough_cnt = (unsigned int***)malloc(sizeof(unsigned int**)*N_AC);
//**  if(hough_cnt==NULL) exit(1);
//**  for(i=0; i<N_AC; i++){
//**    hough_cnt[i] = (unsigned int**)malloc(sizeof(unsigned int*)*DIV_HOUGH_X);
//**    if(hough_cnt[i]==NULL) exit(1);
//**    for(j=0; j<DIV_HOUGH_X; j++){
//**      hough_cnt[i][j] = (unsigned int*)malloc(sizeof(unsigned int)*DIV_HOUGH_Y);
//**      if(hough_cnt[i][j]==NULL) exit(1);
//**    }
//**  }
//**  clear_hough_cnt(hough_cnt);


  //////////////////////////////////////////////////////////////////////////
  neve=0;
  //**  for(theta3_deg=theta_start; theta3_deg<theta_stop; theta3_deg+=THETA3_STEP){
  //  for(itheta3=0; itheta3<N_THETA3; itheta3++){
  //    theta3_deg = rndm->Uniform(THETA3_START, THETA3_STOP);
  for(int ii=0;ii<N_EVENT;){
    printf("%d ",ii);
    theta3_deg = rndm->Uniform(70, 110);
    double Ex2 = 0;//MeV

    if(a.calc(beam_ene,theta3_deg,0,Ex2)==1){
      theta3_deg = rndm->Uniform(70,a.getthr3_max());
      a.calc(beam_ene,theta3_deg,0,Ex2);
    }
    theta3_rad = theta3_deg*(TMath::DegToRad());
    printf("Ex2=%f, theta3=%.2f\n", Ex2, theta3_deg);
    e3 = a.getE_L(2);

//**    for(e3=E3_START; e3<=E3_STOP+0.00001; e3+=E3_STEP){
//**      e3 = 0.7;

    //    for(ie3=0; ie3<N_E3; ie3++){
    //      e3 = rndm->Uniform(E3_START, E3_STOP);
    range = gr_range->Eval(e3);
    range_rec = range + rndm->Gaus(0, RANGE_RESO);
    e3_rec = gr_ene->Eval(range_rec);
    total_cnt=0;
    stop_cnt=0;
    
    srim_alpha->SetKineticEnergy((e3+0.00001)*1.0e6);
    
    /*calculate the kinematics */
    //      ex = calc_Ex(Mass_10C, Mass_4He, Mass_4He, Mass_10C, beam_ene,
    //      		   e3, theta3_rad);
    ex = calc_Ex(Mass_10C, Mass_4He, Mass_4He, Mass_10C, beam_ene,
		 e3_rec, theta3_rad);
    
    theta_cm = calccmang(Mass_10C, Mass_4He, Mass_4He, Mass_10C, beam_ene,
			 ex, theta3_rad, e3_rec);
    
//**      for(phi3_deg=PHI3_START; phi3_deg<=PHI3_STOP+0.00001; phi3_deg+=PHI3_STEP){
    
    phi3_deg = rndm->Uniform(0.,359.);
    phi3_rad = phi3_deg*(TMath::DegToRad());
    dx = range*sin(theta3_rad)*cos(phi3_rad);
    dy = range*sin(theta3_rad)*sin(phi3_rad);	  
    dz = -1.0*range*cos(theta3_rad);
    dr = sqrt(dx*dx+dy*dy+dz*dz);
    
    for(itry=0; itry<N_TRY; itry++){
      neve++;
      total_cnt++;
      stop_inside=0;
      hough_scatt_flag=0;
      
      vtx_x = rndm->Gaus(VTX_X_MEAN, VTX_X_SIGMA);
      vtx_y = rndm->Gaus(VTX_Y_MEAN, VTX_Y_SIGMA);	  
      vtx_z = rndm->Uniform(VTX_Z_START, VTX_Z_STOP);
      
      beam_start_x = vtx_x;
      beam_start_y = vtx_y;
      beam_start_z = 105;	  
      
      stop_x = vtx_x + dx;
      stop_y = vtx_y + dy;
      stop_z = vtx_z + dz;
      
      recoil_track_a_x[0] = vtx_z*N_STRP/TPC_SIZE;
      recoil_track_a_y[0] = vtx_y/driftv;	  
      recoil_track_c_x[0] = vtx_x*N_STRP/TPC_SIZE;
      recoil_track_c_y[0] = vtx_y/driftv;	  
      
      recoil_track_a_x[1] = stop_z*N_STRP/TPC_SIZE;
      recoil_track_a_y[1] = stop_y/driftv;	  
      recoil_track_c_x[1] = stop_x*N_STRP/TPC_SIZE;
      recoil_track_c_y[1] = stop_y/driftv;	  
      
      stop_inside = judge_stop_inside(stop_x, stop_y, stop_z);
      if(range<25) stop_inside = 0;    // range gate added on 2019/01/03
      if(stop_inside==1) stop_cnt++;
      
//**	  h101->Fill(theta3_deg, e3, 1.0);
//**	  if(stop_inside==1) h102->Fill(theta3_deg, e3, 1.0);	  
//**
//**	  h103->Fill(ex, theta_cm);
//**	  if(stop_inside==1) h104->Fill(ex, theta_cm);
//**
//**	  h105->Fill(range, e3);
//**	  if(stop_inside==1) h106->Fill(range, e3);
//**	  
//**	  if(stop_inside==1){
//**	    h200->Fill(vtx_z);
//**	  }
//**
//**	  if(stop_inside==1){
//**	    h201->Fill(range);
//**	    h202->Fill(phi3_deg);
//**	    h203->Fill(phi3_deg, stop_y-vtx_y);
//**	    h204->Fill(e3);
//**	  }

	  /* inject alpha particle for SRIM calc */
	if(stop_inside==1){
	  clear_raw_wave(raw_wave);
	  clear_fadc_data(fadc_data);	    
	  clear_tpc_data(tpc_data);
	  
//**#ifdef VIEW_EVE
//**	    h1000->Reset();
//**	    h1001->Reset();	  
//**	    h1002->Reset();
//**	    h1003->Reset();	  
//**	    h1010->Reset();
//**	    h1011->Reset();	  
//**	    h1012->Reset();
//**	    h1013->Reset();	  
//**	    h1100->Reset();
//**	    h1101->Reset();	    
//**	    h1200->Reset();
//**	    h1201->Reset();	    
//**	    h1300->Reset();
//**	    h1301->Reset();	    
//**	    h2000->Reset();
//**	    h2001->Reset();	    

//**	    for(i=0; i<N_AC; i++){
//**	      for(j=0; j<N_FADC; j++){
//**		hfadc[i][j].Reset();
//**	      }
//**	    }
	    
	  TGraph *gr_alpha_a = new TGraph(2, recoil_track_a_x, recoil_track_a_y);
	  TGraph *gr_alpha_c = new TGraph(2, recoil_track_c_x, recoil_track_c_y);	   
	  
	  gr_alpha_a->SetLineColor(kRed);
	  gr_alpha_c->SetLineColor(kRed);	    
	  gr_alpha_a->SetLineWidth(3);
	  gr_alpha_c->SetLineWidth(3);	    
//**#endif
	  
#ifdef ALPHA_TRACK
	  t_0 = 0.0;
	  if (!srim_alpha->NewTrack(vtx_x*mmTocm, vtx_y*mmTocm, vtx_z*mmTocm, t_0,
				    dx+0.001, dy+0.001, dz+0.001)) {
	    std::cerr << "Generating clusters failed; skipping this track.\n";
	    continue;
	  }	  
	  n_cluster = 0;
	    
	  /* get cluster position of recoil alpha */
	  tot_ne=0;
	  while(srim_alpha->GetCluster(cluster_pos[1], cluster_pos[2], cluster_pos[3],
				       cluster_pos[0],
				       ne, ec, ekin)){
	    
	    n_cluster++;
	    tot_ne+=ne;
	    
//**#ifdef VIEW_EVE
//**	      h1000->Fill(cluster_pos[3]*cmTomm, cluster_pos[2]*cmTomm, ne);
//**	      h1001->Fill(cluster_pos[1]*cmTomm, cluster_pos[2]*cmTomm, ne);
//**	      h1002->Fill(cluster_pos[3]*cmTomm, ne);
//**	      h1003->Fill(cluster_pos[1]*cmTomm, ne);
//**#endif
		
	    if(n_cluster%100==0){
	      printf("cluster:%d, z=%f\n", n_cluster, cluster_pos[3]*cmTomm);
	    }
	    
	    /* drift each electron in the cluster */
#ifdef ONE_ELE
	    int ie=0;
	    int add_ele=1;
#endif

#ifndef ONE_ELE
	    int ie=ne-1;
	    int add_ele=ne;
#endif
	    
	    while(ie<ne){
#ifndef ORIGINAL_DRIFT
	      drift->DriftElectron(cluster_pos[1], cluster_pos[2], cluster_pos[3],
				   cluster_pos[0]);
	      
	      drift->GetElectronEndpoint(0,
					 cluster_pos[1],cluster_pos[2],cluster_pos[3],
					 cluster_pos[0],
					 ele_end_pos[1],ele_end_pos[2],ele_end_pos[3],
					 ele_end_pos[0],
					 drift_status);
#endif
	      
#ifdef ORIGINAL_DRIFT	      
	      drift_electron2(cluster_pos, driftv, diff_tra, diff_long,
			      gen_tra, gen_long, ele_end_pos);
#endif
	      
	      drift_time = (ele_end_pos[0]-cluster_pos[0])*0.10;
//**#ifdef VIEW_EVE
//**		if(ie==0){
//**		  h1010->Fill(ele_end_pos[3]*cmTomm, driftv*drift_time, ne);
//**		  h1011->Fill(ele_end_pos[1]*cmTomm, driftv*drift_time, ne);
//**		  h1012->Fill(ele_end_pos[3]*cmTomm, ne);
//**		  h1013->Fill(ele_end_pos[1]*cmTomm, ne);
//**		}
//**#endif
		
	      //		add_raw_wave(ele_end_pos, drift_time*10.0, add_ele,
	      //			     wave_temp, GAS_GAIN, raw_wave);
	      add_raw_wave2(ele_end_pos, drift_time*10.0, add_ele,
			    wave_spline, GAS_GAIN, raw_wave);
	      ie++;
	    } // end of while(ie<ne)
	  } // end of alpha cluster loop
	  
	  printf("theta=%.1f, e3=%.2f MeV, phi=%.1f, range=%.1f, cluster num=%d, tot_ne=%d\n",
		 theta3_deg, e3, phi3_deg, range, n_cluster, tot_ne);
	  printf("track end position: %.2f, %.2f, %.2f\n",
		 stop_x, stop_y, stop_z);
	  // end of alpha particle track
#endif // end of #ifdef ALPHA_TRACK
	  

//	    /* original clustering */
//	    double tmp_pos[4];
//	    double tmp_ene=e3;
//	    double tmp_stop[4];
//	    double eneloss;
//	    int istep=0;
//	    int nele;
//	    int nele_tot=0;
//	    
//	    tmp_pos[0] = 0;
//	    tmp_pos[1] = vtx_x;
//	    tmp_pos[2] = vtx_y;
//	    tmp_pos[3] = vtx_z;
//
//	    printf("start of original ionization, ene=%f\n", tmp_ene);
//	    while(tmp_ene>W_Val*eVToMeV){
//	      tmp_pos[1]+=dx/dr*ion_step;
//	      tmp_pos[2]+=dy/dr*ion_step;
//	      tmp_pos[3]+=dz/dr*ion_step;	      
//	      eneloss = (gr_eneloss_alpha->Eval(tmp_ene))*eneloss_fac_alpha*ion_step;
//	      tmp_ene-=eneloss;
//	      nele=eneloss/(W_Val*eVToMeV);
//	      nele_tot+=nele;
//	      istep++;
//
//	      /* drift electron*/
//	      for(i=0; i<nele; i++){
//		//		drift_electron2(tmp_pos, driftv, diff_tra, diff_long,
//		//				gen_tra, gen_long, tmp_stop);
//	      }
//	    }
//	    
//	    printf("end of original ionization, ene=%f\n", tmp_ene);
//	    printf("total electron=%d\n", nele_tot);
//	    printf("true stop:(%f,%f,%f), stop:(%f,%f,%f)\n",
//		   stop_x, stop_y, stop_z, tmp_pos[1], tmp_pos[2], tmp_pos[3]);
	    
	    // beam track
#ifdef BEAM_TRACK
	  t_0 = 0.0;
	  if (!srim_10c->NewTrack(beam_start_x*mmTocm, beam_start_y*mmTocm,
				  beam_start_z*mmTocm, t_0,
				  beam_dx+0.0001, beam_dy+0.0001, beam_dz+0.0001)) {
	    std::cerr << "Generating clusters failed; skipping this track.\n";
	    continue;
	  }	  
	  n_cluster = 0;
	  
	  /* get cluster position */
	  while(srim_10c->GetCluster(cluster_pos[1], cluster_pos[2], cluster_pos[3],
				     cluster_pos[0],
				     ne, ec, ekin)){
	    if(n_cluster%1000==0){
	      printf("cluster:%d, z=%f\n", n_cluster, cluster_pos[3]*cmTomm);
	    }
	    if(cluster_pos[3]*cmTomm < -5.0){
	      break;
	    }
	    n_cluster++;

//**#ifdef VIEW_EVE
//**	      h1000->Fill(cluster_pos[3]*cmTomm, cluster_pos[2]*cmTomm, ne);
//**	      h1001->Fill(cluster_pos[1]*cmTomm, cluster_pos[2]*cmTomm, ne);
//**	      h1002->Fill(cluster_pos[3]*cmTomm, ne);
//**	      h1003->Fill(cluster_pos[1]*cmTomm, ne);
//**#endif
		
	    
	    /* drift each electron in the cluster */
#ifdef ONE_ELE
	    int ie=0;
	    int add_ele=1;
#endif
	    
#ifndef ONE_ELE
	    int ie=ne-1;
	    int add_ele=ne;
#endif
	    while(ie<ne){
	      
#ifndef ORIGINAL_DRIFT
	      drift->DriftElectron(cluster_pos[1], cluster_pos[2], cluster_pos[3],
				   cluster_pos[0]);
	      drift->GetElectronEndpoint(0,
					 cluster_pos[1],cluster_pos[2],cluster_pos[3],
					 cluster_pos[0],
					 ele_end_pos[1],ele_end_pos[2],ele_end_pos[3],
					 ele_end_pos[0],
					 drift_status);
#endif
	      
#ifdef ORIGINAL_DRIFT
	      drift_electron2(cluster_pos, driftv, diff_tra, diff_long, gen_tra, gen_long,
			      ele_end_pos);
#endif
	      drift_time = (ele_end_pos[0]-cluster_pos[0])*0.10; 
	      
//**#ifdef VIEW_EVE
//**		if(ie==0){
//**		  h1010->Fill(ele_end_pos[3]*cmTomm, driftv*drift_time, ne);
//**		  h1011->Fill(ele_end_pos[1]*cmTomm, driftv*drift_time, ne);
//**		  h1012->Fill(ele_end_pos[3]*cmTomm, ne);
//**		  h1013->Fill(ele_end_pos[1]*cmTomm, ne);
//**		}
//**#endif
	      
	      //		add_raw_wave(ele_end_pos, drift_time*10.0, add_ele,
	      //			     wave_temp, GAS_GAIN, raw_wave);
	      add_raw_wave2(ele_end_pos, drift_time*10.0, add_ele,
			    wave_spline, GAS_GAIN, raw_wave);
	      ie++;
	    } // end of while(ie<ne)
	  } // end of 10c cluster loop
	  
#endif  // end of #ifdev BEAM_TRACK
	  make_fadc_data(raw_wave, fadc_data);
	  add_fadc_noise(fadc_data, gen_fadc_noise, FADC_NOISE);
	  make_tpc_data(raw_wave, tpc_data, TPC_THRESHOLD);

	  int filenum = ii%filesize;
	  for(int jj=0;jj<N_AC;jj++){
	    for(int kk=0;kk<N_TCLK;kk++){
	      for(int ll=0;ll<N_STRP;ll++){
		fprintf(fp_t[filenum],"%d ",tpc_data[jj][kk][ll]);
	      }
	    }
	  }
	  fprintf(fp_t[filenum],"\n");
	  fprintf(fp_r[filenum],"%f %f %f %f %f %d %d %d %d %d %d %d %d\n",
		  theta3_deg,phi3_deg,range,e3,Ex2,
		  int(recoil_track_a_x[0]),int(recoil_track_a_y[0]),int(recoil_track_a_x[1]),int(recoil_track_a_y[1]),
		  int(recoil_track_c_x[0]),int(recoil_track_c_y[0]),int(recoil_track_c_x[1]),int(recoil_track_c_y[1]));
	  ii++;
	  
//**#ifdef VIEW_EVE
//**	    for(j=0; j<N_TCLK; j++){
//**	      for(k=0; k<N_STRP; k++){
//**		h1100->SetBinContent(k,j, raw_wave[0][j][k]);
//**		h1101->SetBinContent(k,j, raw_wave[1][j][k]);		  
//**		if(tpc_data[0][j][k]>0) h1300->SetBinContent(k,j, tpc_data[0][j][k]);
//**		if(tpc_data[1][j][k]>0) h1301->SetBinContent(k,j, tpc_data[1][j][k]);		  
//**	      }
//**	    }
//**	    for(j=0; j<N_FADC; j++){
//**	      for(k=0; k<N_ACLK; k++){
//**		h1200->Fill(j,k, fadc_data[0][j][k]);
//**		h1201->Fill(j,k, fadc_data[1][j][k]);		  
//**	      }
//**	    }
//**
//**	    for(i=0; i<N_AC; i++){
//**	      for(j=0; j<N_FADC; j++){
//**		for(k=0; k<N_ACLK; k++){
//**		  hfadc[i][j].SetBinContent(k,fadc_data[i][j][k]);
//**		}
//**	      }
//**	    }
//**	    
//**	    c1->cd();
//**	    c1->Clear();
//**	    c1->Divide(2,1);
//**	    c1->cd(1);
//**	    h1000->Draw("colz");
//**	    c1->cd(2);
//**	    h1001->Draw("colz");
//**	    c1->Update();
//**
//**	    c2->cd();
//**	    c2->Clear();
//**	    c2->Divide(2,1);
//**	    c2->cd(1);
//**	    h1002->Draw("");
//**	    c2->cd(2);
//**	    h1003->Draw("");
//**	    c2->Update();
//**
//**	    c3->cd();
//**	    c3->Clear();
//**	    c3->Divide(2,1);
//**	    c3->cd(1);
//**	    h1010->Draw("colz");
//**	    c3->cd(2);
//**	    h1011->Draw("colz");
//**	    c3->Update();
//**
//**	    c4->cd();
//**	    c4->Clear();
//**	    c4->Divide(2,1);
//**	    c4->cd(1);
//**	    h1012->Draw("");
//**	    c4->cd(2);
//**	    h1013->Draw("");
//**	    c4->Update();
//**
//**	    c5->cd();
//**	    c5->Clear();
//**	    c5->Divide(2,1);
//**	    c5->cd(1);
//**	    h1100->Draw("colz");
//**	    c5->cd(2);
//**	    h1101->Draw("colz");
//**	    c5->Update();
//**
//**	    c6->cd();
//**	    c6->Clear();
//**	    c6->Divide(2,1);
//**	    c6->cd(1);
//**	    h1300->Draw("box");
//**	    gr_alpha_a->Draw("L");
//**	    c6->cd(2);
//**	    h1301->Draw("box");
//**	    gr_alpha_c->Draw("L");
//**	    c6->Update();
//**#endif
//**
//**	    ///////////////////////////
//**	    //   tracking analysis   //
//**	    ///////////////////////////
//**	    printf("\n");
//**	    printf("******* Start of Analysis ***********\n");
//**	    
//**	    /* initialization of parameters */
//**	    clear_hough_cnt(hough_cnt);
//**	    for(i=0; i<N_AC; i++){
//**	      hough_max_pos[i][0]=0;
//**	      hough_max_pos[i][1]=0;	      
//**	      hough_max_cnt[i]=0;
//**	    }
//**
//**	    for(i=0; i<N_AC; i++){
//**	      for(j=0; j<N_FADC; j++){
//**		for(k=0; k<FADC_MAX_PULSE; k++){
//**		  pulse_integ[i][j][k]=0;
//**		  pulse_lead[i][j][k]=0;
//**		  pulse_width[i][j][k]=0;
//**		}
//**	      }
//**	    }
//**	    /* initialization to here */	    
//**	    
//**
//**	    /* get FADC pulse */
//**	    calc_pulse_integ_fadc(fadc_data, pulse_integ, pulse_lead,
//**				  pulse_width, FADC_PULSE_TH);
//**
//**	    for(i=0; i<N_FADC; i++){
//**	      for(j=0; j<FADC_MAX_PULSE; j++){
//**		if(pulse_integ[ANODE][i][j]>0){
//**		  printf("ch=%d, index=%d, integ=%d, lead=%d, width=%d\n",
//**			 i,j, pulse_integ[ANODE][i][j], pulse_lead[ANODE][i][j],
//**			 pulse_width[ANODE][i][j]);
//**		}
//**	      }
//**	    }
//**	    
//**	    
//**	    /* transform into Hough space */
//**	    hough_tra(tpc_data, hough_cnt);
//**	    
//**	    /* find max position */
//**	    find_hough_max(hough_cnt, hough_max_pos, hough_max_cnt);
//**	    for(i=0; i<N_AC; i++){
//**	      hough_theta[i] = (hough_max_pos[i][0]*360.0)/DIV_HOUGH_X;
//**	      hough_r[i] = (hough_max_pos[i][1]*1024.0)/DIV_HOUGH_Y;
//**	    }
//**
//**	    printf("--------------------------------------\n");
//**	    printf("Anode Hough: pos:(%.1f, %.1f), cnt:%d\n",
//**		   hough_theta[0], hough_r[0], hough_max_cnt[0]);
//**
//**	    printf("Cathode Hough: pos:(%.1f, %.1f), cnt:%d\n",
//**		   hough_theta[1], hough_r[1], hough_max_cnt[1]);
//**		   
//**	    for(i=0; i<N_AC; i++){
//**	      line_para[i][0]=hough_r[i]/sin(hough_theta[i]*TMath::DegToRad());
//**	      line_para[i][1]=-1.0/tan(hough_theta[i]*TMath::DegToRad());
//**
//**	      line_y[i][0]=line_para[i][0];
//**	      line_y[i][1]=line_para[i][0]+N_STRP*line_para[i][1];
//**	    }
//**	    
//**
//**	    /* scattering detection */
//**	    if(hough_max_cnt[ANODE]   > HOUGH_TH_A &&
//**	       hough_max_cnt[CATHODE] > HOUGH_TH_C &&
//**	       (hough_theta[ANODE]<(90-HOUGH_RECOIL_ANG_A) ||
//**		hough_theta[ANODE]>(90+HOUGH_RECOIL_ANG_A)) ){
//**	      hough_scatt_flag=1;
//**	      printf("scattering!\n");
//**	    }
//**	    
//**	    if(stop_inside==1 && hough_scatt_flag==1){
//**	      h300->Fill(theta3_deg, e3, 1.0);
//**	      h301->Fill(ex, theta_cm, 1.0);
//**	      h302->Fill(phi3_deg, 1.0);
//**	      h303->Fill(vtx_z, 1.0);
//**	    }

//**#ifdef VIEW_EVE
//**	    /* fill into histograms */
//**	    for(i=0; i<N_AC; i++){
//**	      for(j=0; j<DIV_HOUGH_X; j++){
//**		for(k=0; k<DIV_HOUGH_Y; k++){
//**		  if(i==0){
//**		    h2000->Fill(j*360.0/DIV_HOUGH_X, k*1024.0/DIV_HOUGH_Y,
//**				hough_cnt[i][j][k]);
//**		  }
//**		  if(i==1){
//**		    h2001->Fill(j*360.0/DIV_HOUGH_X, k*1024.0/DIV_HOUGH_Y,
//**				hough_cnt[i][j][k]);
//**		  }
//**
//**		}
//**	      }
//**	    }
//**
//**	    TGraph *gr_hough_a = new TGraph();
//**	    gr_hough_a->SetPoint(0,   0.0,line_y[0][0]);
//**	    gr_hough_a->SetPoint(1, 256.0,line_y[0][1]);	    
//**	    gr_hough_a->SetLineColor(3);
//**	    gr_hough_a->SetLineWidth(2);
//**	    
//**	    TGraph *gr_hough_c = new TGraph();	    
//**	    gr_hough_c->SetPoint(0,   0.0,line_y[1][0]);
//**	    gr_hough_c->SetPoint(1, 256.0,line_y[1][1]);	    
//**	    gr_hough_c->SetLineColor(3);
//**	    gr_hough_c->SetLineWidth(2);
//**
//**	    c10->cd();
//**	    c10->Clear();
//**	    c10->Divide(2,1);
//**	    c10->cd(1);
//**	    h2000->Draw("colz");
//**	    c10->cd(2);
//**	    h2001->Draw("colz");
//**	    c10->Update();
//**
//**	    c6->cd(1);
//**	    gr_hough_a->Draw("L");
//**	    c6->cd(2);
//**	    gr_hough_c->Draw("L");
//**	    c6->Update();
//**	    
//**	    TFile *testfile = new TFile("tmp.root", "recreate");
//**	    h1100->Write();
//**	    h1101->Write();
//**	    h1300->Write();
//**	    h1301->Write();
//**	    h2000->Write();
//**	    h2001->Write();	    
//**	    for(i=0; i<N_AC; i++){
//**	      for(j=0; j<N_FADC; j++){
//**		hfadc[i][j].Write();
//**	      }
//**	    }
//**	    testfile->Close();


//**	  while(getchar() != '\n');
//**#endif
	  
	} // end of if(stop_inside==1)
	
      }  // end of itry loop
//**      }  // end of phi3 loop
      
      //      printf("theta=%.2f, e3=%.2f, eff=%.4f\n", theta3_deg, e3, stop_eff);
      stop_eff = stop_cnt/((double)total_cnt);
//**      h100->Fill(theta3_deg, e3, stop_eff);
      
//**    }  // end of e3 loop
  }// end of ii loop
//**  }  // end of theta3 loop

  /* save dat file */
  for(int ii=0;ii<filesize;ii++){
    fclose(fp_t[ii]);
    fclose(fp_r[ii]);
  }

//**  /* output root file */
//**  TFile *outfile = new TFile(outfname, "recreate");
//**  h100->Write();
//**  h101->Write();
//**  h102->Write();
//**  h103->Write();
//**  h104->Write();
//**  h105->Write();
//**  h106->Write();    
//**  h200->Write();
//**  h201->Write();
//**  h202->Write();
//**  h203->Write();
//**  h204->Write();    
//**  h300->Write();
//**  h301->Write();
//**  h302->Write();
//**  h303->Write();  
//**  outfile->Close();
//**
//**  printf("total event: %d\n", neve);

//**#ifdef VIEW_EVE
//**  app.Run();
//**#endif

  ///////////////////////////////
  //   release memory
  ///////////////////////////////
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      free(raw_wave[i][j]);
    }
    free(raw_wave[i]);
  }
  free(raw_wave);

  for(i=0; i<N_AC; i++){
    for(j=0; j<N_FADC; j++){
      free(fadc_data[i][j]);
    }
    free(fadc_data[i]);
  }
  free(fadc_data);
  
  for(i=0; i<N_AC; i++){
    for(j=0; j<N_TCLK; j++){
      free(tpc_data[i][j]);
    }
    free(tpc_data[i]);
  }
  free(tpc_data);

//**  for(i=0; i<N_AC; i++){
//**    for(j=0; j<DIV_HOUGH_X; j++){
//**      free(hough_cnt[i][j]);
//**    }
//**    free(hough_cnt[i]);
//**  }
//**  free(hough_cnt);

  return 0;
}


/* judge whether the recoil stoped inside the TPC volume */
int judge_stop_inside(double stop_x, double stop_y, double stop_z){
  int stop=0;
  int good_x=0;
  int good_y=0;
  int good_z=0;  

  if((stop_x>STOP_X_MIN1 && stop_x<STOP_X_MAX1) ||
     (stop_x>STOP_X_MIN2 && stop_x<STOP_X_MAX2)){
    good_x=1;
  }
  if(stop_y>STOP_Y_MIN && stop_y<STOP_Y_MAX) good_y=1;
  if(stop_z>STOP_Z_MIN && stop_z<STOP_Z_MAX) good_z=1;  

  if(good_x==1 && good_y==1 && good_z==1) stop=1;
  
  return stop;
}

