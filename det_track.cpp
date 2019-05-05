#define SRIM_PARA
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
#include <TArrow.h>

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

#include "para.h"
#include "anatpc.h"
#include "kinema.h"
#include "mkdata.h"
#include "12C_3alpha.hpp"

#define VIEW_EVE 1
#define ALPHA_TRACK 1
#define SAVE_FILE 1
#ifndef SAVE_FILE
//#define DRAW_ARROW 1
//#define DRAW_CANVAS 1
#endif
//#define BEAM_TRACK 1
//#define ORIGINAL_DRIFT 1
//#define ONE_ELE 1

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
  
    
  /* initialize random seed (always same init) */
  TRandom3 *rndm = new TRandom3();
  rndm->SetSeed(index);
  
  /* read energy to range table */
  char grfilename[512], gr2filename[512];
  sprintf(grfilename,  "%s/tables/ene_to_range_%s_%d.dat", workdir, detection_gas, press);
  sprintf(gr2filename, "%s/tables/range_to_ene_%s_%d.dat", workdir, detection_gas, press);  
  printf("range table: %s\n", grfilename);
  TGraph *gr_range;
  gr_range = new TGraph(grfilename);
  TGraph *gr_ene;
  gr_ene = new TGraph(gr2filename);
  
  /* read energy loss table */
//****  char enelossfname_alpha[512];
//****  sprintf(enelossfname_alpha, "%s/tables/eneloss_alpha_%d.dat", workdir, press);
//****  TGraph *gr_eneloss_alpha = new TGraph(enelossfname_alpha);
  
  TLorentzVector alpha[nAlpha];
  double phi3_deg[nAlpha];
  double phi3_rad[nAlpha];     // recoil phi angle (radian)
  double theta3_deg[nAlpha];
  double theta3_rad[nAlpha];   // recoil angle in LAB (radian)  
  double e3[nAlpha];           // recoil energy in LAB (MeV)
  double e3_rec[nAlpha];       // reconstructed recoil energy with finite resolution in LAB (MeV)  
  double range[nAlpha];        // recoil range (mm)
  double range_rec[nAlpha];    // reconstructed recoil range with finite resolution (mm)  
  double vtx_x,vtx_y,vtx_z;     // vertex position (mm)
  double dx[nAlpha],dy[nAlpha],dz[nAlpha];              // position difference between vertex and stop points
  double dr[nAlpha];
  double stop_x[nAlpha],stop_y[nAlpha],stop_z[nAlpha];  // recoil particle stop position (mm)
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

  double recoil_track_a_x[nAlpha][2], recoil_track_a_y[nAlpha][2];
  double recoil_track_c_x[nAlpha][2], recoil_track_c_y[nAlpha][2];  
  
  /* Hough analysis parameters */
  int pulse_integ[N_AC][N_FADC][FADC_MAX_PULSE];
  int pulse_lead[N_AC][N_FADC][FADC_MAX_PULSE];
  int pulse_width[N_AC][N_FADC][FADC_MAX_PULSE];  
  
  double line_para[N_AC][2];
  double line_y[N_AC][2];

  /* Event viewer */
  TApplication app("app",&argc,argv);
  TH2D *h_anode = new TH2D("h_anode","h_anode",256,0.,256.,1024,0.,1024.);
  TH2D *h_cathode = new TH2D("h_cathode","h_cathode",256,0.,256.,1024,0.,1024.);
  TH2D *h_raw_anode = new TH2D("h_raw_anode","h_raw_anode",256,0.,256.,1024,0.,1024.);
  TH2D *h_raw_cathode = new TH2D("h_raw_cathode","h_raw_cathode",256,0.,256.,1024,0.,1024.);
#ifdef DRAW_CANVAS
  TCanvas *c1 = new TCanvas("c1");
  TCanvas *c2 = new TCanvas("c2");
  TCanvas *c3 = new TCanvas("c3");
  TCanvas *c4 = new TCanvas("c4");
#endif
#ifdef DRAW_ARROW
  TArrow ar_anode[3];
  TArrow ar_cathode[3];
#endif
#ifdef SAVE_FILE
  std::string tempname = FILENAME+"tpc.dat";
  FILE *file_tpc = fopen(tempname.c_str(),"w");
  tempname = FILENAME+"para.dat";
  FILE *file_para = fopen(tempname.c_str(),"w");
#endif

  
  ///////////////////////////
  //   garfield settings   //
  ///////////////////////////

  /* Read Magboltz file */
  char magfname[512];
  sprintf(magfname, "%s/tables/%s_%d.gas", workdir, detection_gas, press);

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
    gas->GetElectronVelocityE(i, 0, 0, velocity[i]); // cm/ns
    gas->GetElectronTransverseDiffusion(i, 0, 0, diff_t[i]); // sqrt(cm)
    gas->GetElectronLongitudinalDiffusion(i, 0, 0, diff_l[i]); // sqrt(cm)
    gas->GetElectronTownsend(i, 0, 0, town_a[i]);
    gas->GetElectronAttachment(i, 0, 0, attach[i]);
    
    ef[i] = efields[i];
//****    velocity[i] = velocity[i]*1000; // 10um/ns
//****    diff_t[i] = diff_t[i]*10000; // sqrt(mm)
//****    diff_l[i] = diff_l[i]*10000; // sqrt(mm)
    
    if(town_a[i]<0) town_a[i] = 0;
    if(attach[i]<0) attach[i] = 0;
  }  

  /* calculate drift velocity */
  TGraph *gr_driftv = new TGraph(efields.size(), ef, velocity);
  double driftv = gr_driftv->Eval(E_FIELD)*10; // mm/ns
  
  /* calculate trans diffusion */ 
  TGraph *gr_diff_tra = new TGraph(efields.size(), ef, diff_t);
  double diff_tra = gr_diff_tra->Eval(E_FIELD)*10000.; //sqrt(mm)

  /* calculate long diffusion */ 
  TGraph *gr_diff_long = new TGraph(efields.size(), ef, diff_l);
  double diff_long = gr_diff_long->Eval(E_FIELD)*10000.; //sqrt(mm)

  printf("E: %.1f V/cm, driftv: %.3f mm/ns, diff_tra: %.1f, long_diff: %.1f\n",
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
//****  char srimfname_alpha[512], srimfname_10c[512];
  char srimfname_alpha[512];
  sprintf(srimfname_alpha, "%s/tables/4He_%s_%d.srim", workdir, detection_gas, press);
//****  sprintf(srimfname_10c,   "%s/tables/10C_HeCO2_96_4_%d.srim", workdir, press);  

  TrackSrim* srim_alpha = new TrackSrim();
  srim_alpha->SetSensor(sensor);
  if(srim_alpha->ReadFile(srimfname_alpha)==0){
    printf("error in reading SRIM file: %s\n", srimfname_alpha);
    exit(0);
  }

  srim_alpha->SetWorkFunction(W_Val[detection_gas]);  // in eV
  srim_alpha->SetFanoFactor(Fano_Factor);
  srim_alpha->SetModel(4);
  srim_alpha->SetAtomicMassNumbers(Mass[detection_gas],Charge[detection_gas]); // changed by T. Doi 2019/3/28
  srim_alpha->SetDensity(density[detection_gas]);

  srim_alpha->SetTargetClusterSize(Cluster_Size);  
  srim_alpha->DisableTransverseStraggling();
  srim_alpha->DisableLongitudinalStraggling();

  printf("************************************************************\n");
  srim_alpha->Print();
  printf("************************************************************\n");

  // random generators
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

#ifdef SAVE_FILE
//****  char branch[1024];
//****  tree->Branch("h_anode","TH2D",h_anode,128000,0);
//****  tree->Branch("h_cathode","TH2D",h_cathode,128000,0);
//****  tree->Branch("h_raw_anode","TH2D",h_raw_anode,128000,0);
//****  tree->Branch("h_raw_cathode","TH2D",h_raw_cathode,128000,0);
//****  tree->Branch("raw_wave",raw_wave,"raw_wave/D");
//****  tree->Branch("fadc_data",fadc_data,"fadc_data/I");
//****  tree->Branch("tpc_data",tpc_data,"tpc_data/I");
//****  sprintf(branch,"recoil_track_a_x[%d][2]/D",nAlpha);
//****  tree->Branch("recoil_track_a_x",recoil_track_a_x,branch);
//****  sprintf(branch,"recoil_track_a_y[%d][2]/D",nAlpha);
//****  tree->Branch("recoil_track_a_y",recoil_track_a_y,branch);
//****  sprintf(branch,"recoil_track_c_x[%d][2]/D",nAlpha);
//****  tree->Branch("recoil_track_c_x",recoil_track_c_x,branch);
//****  sprintf(branch,"recoil_track_c_y[%d][2]/D",nAlpha);
//****  tree->Branch("recoil_track_c_y",recoil_track_c_y,branch);
//****  sprintf(branch,"theta3_deg[%d]/D",nAlpha);
//****  tree->Branch("theta3_deg",theta3_deg,branch);
//****  sprintf(branch,"phi3_deg[%d]/D",nAlpha);
//****  tree->Branch("phi3_deg",phi3_deg,branch);
//****  sprintf(branch,"e3[%d]/D",nAlpha);
//****  tree->Branch("e3",e3,branch);
//****  sprintf(branch,"dr[%d]/D",nAlpha);
//****  tree->Branch("dr",dr,branch);
#endif
  
  
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

//****  const int filesize = N_EVENT/MAX_FILE_EVENT;
//****  FILE *fp_t[filesize];
//****  FILE *fp_r[filesize];
//****  char filename[1024];
//****  for(int ii=0;ii<filesize;ii++){
//****    sprintf(filename,"simu_track-%d.dat",ii);
//****    fp_t[ii] = fopen(filename,"w");
//****    sprintf(filename,"simu_result-%d.dat",ii);
//****    fp_r[ii] = fopen(filename,"w");
//****  }


  //////////////////////////////////////////////////////////////////////////
  neve=0;
  TStyle *style = new TStyle("Plain","style1");
//****  style->SetOptStat("ourme");
  style->cd();
  //**  for(theta3_deg=theta_start; theta3_deg<theta_stop; theta3_deg+=THETA3_STEP){
  //  for(itheta3=0; itheta3<N_THETA3; itheta3++){
  //    theta3_deg = rndm->Uniform(THETA3_START, THETA3_STOP);
  for(int ii=0;ii<N_EVENT;){
    printf("%d \n",ii);

    // set particle data
    decay12C(alpha);
//    theta3_deg = alpha[0].Theta();
//    theta3_rad = theta3_deg*(TMath::DegToRad());
//    phi3_deg = rndm->Uniform(0.,359.);
//    phi3_rad = phi3_deg*(TMath::DegToRad());
#ifdef ALPHA_TRACK
    
    neve++;
    total_cnt++;
    stop_inside=0;
    
    vtx_x = rndm->Gaus(VTX_X_MEAN, VTX_X_SIGMA);
    vtx_y = rndm->Gaus(VTX_Y_MEAN, VTX_Y_SIGMA);  
    vtx_z = rndm->Uniform(VTX_Z_START, VTX_Z_STOP);
    
    int alpha_inside = 0;
    
    clear_raw_wave(raw_wave);
    clear_fadc_data(fadc_data);	    
    clear_tpc_data(tpc_data);
    
    for(int i_alpha=0;i_alpha<nAlpha;i_alpha++){
      theta3_rad[i_alpha] = alpha[i_alpha].Theta();
      phi3_rad[i_alpha] = alpha[i_alpha].Phi();
      theta3_deg[i_alpha] = theta3_rad[i_alpha]*(TMath::RadToDeg());
      phi3_deg[i_alpha] = phi3_rad[i_alpha]*(TMath::RadToDeg());
      //    printf("Ex2=%f, theta3=%.2f\n", Ex2, theta3_deg);
      e3[i_alpha] = (alpha[i_alpha].E()-M_alpha)*MeV;
      
      range[i_alpha] = gr_range->Eval(e3[i_alpha]);
      range_rec[i_alpha] = range[i_alpha] + rndm->Gaus(0, RANGE_RESO);
      e3_rec[i_alpha] = gr_ene->Eval(range_rec[i_alpha]);
      total_cnt=0;
      stop_cnt=0;
      
//****      srim_alpha->SetKineticEnergy((e3+0.00001)*1.0e6);
      
      /*calculate the kinematics */
      //      ex = calc_Ex(Mass_10C, Mass_4He, Mass_4He, Mass_10C, beam_ene,
      //      		   e3, theta3_rad);
//****      ex = calc_Ex(Mass_10C, Mass_4He, Mass_4He, Mass_10C, beam_ene,
//****		   e3_rec, theta3_rad);
//****      
//****      theta_cm = calccmang(Mass_10C, Mass_4He, Mass_4He, Mass_10C, beam_ene,
//****			   ex, theta3_rad, e3_rec);
      
//****      dx = range*sin(theta3_rad)*cos(phi3_rad);
//****      dy = range*sin(theta3_rad)*sin(phi3_rad);	  
//****      dz = -1.0*range*cos(theta3_rad);
//****      dr = sqrt(dx*dx+dy*dy+dz*dz);
      dx[i_alpha] = range[i_alpha]*alpha[i_alpha].X()/sqrt(alpha[i_alpha].X()*alpha[i_alpha].X()+
							   alpha[i_alpha].Y()*alpha[i_alpha].Y()+
							   alpha[i_alpha].Z()*alpha[i_alpha].Z());
      dy[i_alpha] = range[i_alpha]*alpha[i_alpha].Y()/sqrt(alpha[i_alpha].X()*alpha[i_alpha].X()+
							   alpha[i_alpha].Y()*alpha[i_alpha].Y()+
							   alpha[i_alpha].Z()*alpha[i_alpha].Z());
      dz[i_alpha] = -1.0*range[i_alpha]*alpha[i_alpha].Z()/sqrt(alpha[i_alpha].X()*alpha[i_alpha].X()+
								alpha[i_alpha].Y()*alpha[i_alpha].Y()+
								alpha[i_alpha].Z()*alpha[i_alpha].Z());
      dr[i_alpha] = sqrt(dx[i_alpha]*dx[i_alpha]+dy[i_alpha]*dy[i_alpha]+dz[i_alpha]*dz[i_alpha]);
      
      stop_x[i_alpha] = vtx_x + dx[i_alpha];
      stop_y[i_alpha] = vtx_y + dy[i_alpha];
      stop_z[i_alpha] = vtx_z + dz[i_alpha];
    } // end of for(int i_alpha=...)

    double trig_time = GetTrigTime(vtx_y, stop_y, driftv);
    
    for(int i_alpha=0;i_alpha<nAlpha;i_alpha++){  
      
//****      recoil_track_a_x[i_alpha][0] = vtx_z*N_STRP/TPC_SIZE;
//****      recoil_track_a_y[i_alpha][0] = (vtx_y/driftv+trig_time)*SAMPLING_RATIO; // mm*ns/mm = ns/freq(GHz) = ch
//****      recoil_track_c_x[i_alpha][0] = vtx_x*N_STRP/TPC_SIZE;
//****      recoil_track_c_y[i_alpha][0] = (vtx_y/driftv+trig_time)*SAMPLING_RATIO;	  
//****      
//****      recoil_track_a_x[i_alpha][1] = stop_z[i_alpha]*N_STRP/TPC_SIZE;
//****      recoil_track_a_y[i_alpha][1] = (stop_y[i_alpha]/driftv+trig_time)*SAMPLING_RATIO;
//****      recoil_track_c_x[i_alpha][1] = stop_x[i_alpha]*N_STRP/TPC_SIZE;
//****      recoil_track_c_y[i_alpha][1] = (stop_y[i_alpha]/driftv+trig_time)*SAMPLING_RATIO;
//****      ar_anode[i_alpha] = TArrow(recoil_track_a_x[i_alpha][0],recoil_track_a_y[i_alpha][0],
//****				 recoil_track_a_x[i_alpha][1],recoil_track_a_y[i_alpha][1],0.01,"|>");
//****      ar_cathode[i_alpha] = TArrow(recoil_track_c_x[i_alpha][0],recoil_track_c_y[i_alpha][0],
//****				   recoil_track_c_x[i_alpha][1],recoil_track_c_y[i_alpha][1],0.01,"|>");
//****      ar_anode[i_alpha].SetLineColor(kRed);
//****      ar_cathode[i_alpha].SetLineColor(kRed);
//****      ar_anode[i_alpha].SetFillColor(kRed);
//****      ar_cathode[i_alpha].SetFillColor(kRed);
//****      ar_anode[i_alpha].SetLineWidth(3);
//****      ar_cathode[i_alpha].SetLineWidth(3);

      stop_inside = judge_stop_inside(stop_x[i_alpha], stop_y[i_alpha], stop_z[i_alpha]);
//****      if(range<25) stop_inside = 0;    // range gate added on 2019/01/03
//****      if(stop_inside==1){
//****
//****      }else{
//****	std::cout << "alpha particle is out of range" << std::endl;
//****	std::cout << "E_alpha = " << e3 << " ";
//****	std::cout << "stop_x = " << stop_x << " ";
//****	std::cout << "stop_y = " << stop_y << " ";
//****	std::cout << "stop_z = " << stop_z << std::endl;
//****      }
      
      /* inject alpha particle for SRIM calc */
      srim_alpha->SetKineticEnergy(e3[i_alpha]*1.0e6);
      if(stop_inside==1){
	stop_cnt++;
	alpha_inside++;
	
	t_0 = 0.0;
//****	if (!srim_alpha->NewTrack(vtx_x*mmTocm, vtx_y*mmTocm, vtx_z*mmTocm, t_0,
//****				  dx+0.001, dy+0.001, dz+0.001)) {
//****	  std::cerr << "Generating clusters failed; skipping this track.\n";
//****	  continue;
//****	}
	if(!srim_alpha->NewTrack(vtx_x*mmTocm, vtx_y*mmTocm, vtx_z*mmTocm, t_0, dx[i_alpha], dy[i_alpha], dz[i_alpha])){
	  std::cerr << "Generating clusters failed; skipping this track." << std::endl;
	  continue;
	}
	n_cluster = 0;
	
	/* get cluster position of recoil alpha */
	tot_ne=0;
	
	bool first_flag = true;
	double vtx_drift_time;
	
	while(srim_alpha->GetCluster(cluster_pos[1], cluster_pos[2], cluster_pos[3],
				     cluster_pos[0],
				     ne, ec, ekin)){
	  
	  n_cluster++;
	  tot_ne+=ne;
	  
	  
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
	    if(drift->AvalancheElectron(cluster_pos[1], cluster_pos[2], cluster_pos[3], cluster_pos[0])){
	      double ne_sub = drift->GetNumberOfElectronEndpoints();
	      for(int ie_sub=0;ie_sub<ne_sub;ie_sub++){
		drift->GetElectronEndpoint(ie_sub,// original is 0
					   cluster_pos[1],cluster_pos[2],cluster_pos[3],
					   cluster_pos[0],
					   ele_end_pos[1],ele_end_pos[2],ele_end_pos[3],
					   ele_end_pos[0],
					   drift_status);
		drift_time = (ele_end_pos[0]-cluster_pos[0]);
		//****		drift_time = ele_end_pos[0]-trig_time;
		add_raw_wave2(ele_end_pos, drift_time, add_ele,
			      wave_spline, GAS_GAIN, raw_wave);
		if(first_flag){
		  vtx_drift_time = drift_time;
		  first_flag = false;
		}
	      } // end of for(int ie_sub...)
	    }// end of if(drift->...)
	    ie++;
#endif

#ifdef ORIGINAL_DRIFT	      
	    drift_electron2(cluster_pos, driftv, diff_tra, diff_long,
			    gen_tra, gen_long, ele_end_pos);
	    
//****	    drift_time = (ele_end_pos[0]-cluster_pos[0]);
	    drift_time = ele_end_pos[0]-trig_time;
	    
	    //		add_raw_wave(ele_end_pos, drift_time, add_ele,
	    //			     wave_temp, GAS_GAIN, raw_wave);
	    add_raw_wave2(ele_end_pos, drift_time, add_ele,
			  wave_spline, GAS_GAIN, raw_wave);
	    ie++;
#endif
	  } // end of while(ie<ne)
	} // end of alpha cluster loop
	
	recoil_track_a_x[i_alpha][0] = vtx_z*N_STRP/TPC_SIZE;
	recoil_track_a_y[i_alpha][0] = (vtx_drift_time+BUFF_TIME)*SAMPLING_RATIO; // mm*ns/mm = ns/freq(GHz) = ch
	recoil_track_c_x[i_alpha][0] = vtx_x*N_STRP/TPC_SIZE;
	recoil_track_c_y[i_alpha][0] = (vtx_drift_time+BUFF_TIME)*SAMPLING_RATIO;	  
	
	recoil_track_a_x[i_alpha][1] = cluster_pos[3]*cmTomm*N_STRP/TPC_SIZE;
	recoil_track_a_y[i_alpha][1] = (drift_time+BUFF_TIME)*SAMPLING_RATIO;
	recoil_track_c_x[i_alpha][1] = cluster_pos[1]*cmTomm*N_STRP/TPC_SIZE;
	recoil_track_c_y[i_alpha][1] = (drift_time+BUFF_TIME)*SAMPLING_RATIO;
#ifdef DRAW_ARROW
	ar_anode[i_alpha] = TArrow(recoil_track_a_x[i_alpha][0],recoil_track_a_y[i_alpha][0],
				   recoil_track_a_x[i_alpha][1],recoil_track_a_y[i_alpha][1],0.01,"|>");
	ar_cathode[i_alpha] = TArrow(recoil_track_c_x[i_alpha][0],recoil_track_c_y[i_alpha][0],
				     recoil_track_c_x[i_alpha][1],recoil_track_c_y[i_alpha][1],0.01,"|>");
	ar_anode[i_alpha].SetLineColor(kRed);
	ar_cathode[i_alpha].SetLineColor(kRed);
	ar_anode[i_alpha].SetFillColor(kRed);
	ar_cathode[i_alpha].SetFillColor(kRed);
	ar_anode[i_alpha].SetLineWidth(3);
	ar_cathode[i_alpha].SetLineWidth(3);
#endif

	
//****	printf("theta=%.1f, e3=%.2f MeV, phi=%.1f, range=%.1f, cluster num=%d, tot_ne=%d\n",
//****	       theta3_rad, e3, phi3_rad, range, n_cluster, tot_ne);
	printf("theta=%.1f, e3=%.2f MeV, phi=%.1f, range=%.1f, cluster num=%d, tot_ne=%d\n",
	       theta3_deg[i_alpha], e3[i_alpha], phi3_deg[i_alpha], range[i_alpha], n_cluster, tot_ne);
	printf("track end position: %.2f, %.2f, %.2f\n",
	       stop_x[i_alpha], stop_y[i_alpha], stop_z[i_alpha]);
	// end of alpha particle track
	
#endif // end of #ifdef ALPHA_TRACK
	
	
      
//****	int filenum = ii%filesize;
//****	for(int jj=0;jj<N_AC;jj++){
//****	  for(int kk=0;kk<N_TCLK;kk++){
//****	    for(int ll=0;ll<N_STRP;ll++){
//****	      fprintf(fp_t[filenum],"%d ",tpc_data[jj][kk][ll]);
//****	    }
//****	  }
//****	}
//****	fprintf(fp_t[filenum],"\n");
//****	fprintf(fp_r[filenum],"%f %f %f %f %f %d %d %d %d %d %d %d %d\n",
//****		theta3_rad,phi3_rad,range,e3,ex,
//****		int(recoil_track_a_x[0]),int(recoil_track_a_y[0]),int(recoil_track_a_x[1]),int(recoil_track_a_y[1]),
//****		int(recoil_track_c_x[0]),int(recoil_track_c_y[0]),int(recoil_track_c_x[1]),int(recoil_track_c_y[1]));
  
      } // end of if(stop_inside==1)
    }  // end of i_alpha loop

    if(alpha_inside == nAlpha){
      make_fadc_data(raw_wave, fadc_data);
      add_fadc_noise(fadc_data, gen_fadc_noise, FADC_NOISE);
      make_tpc_data(raw_wave, tpc_data, TPC_THRESHOLD);
      for(int clk_id=0;clk_id<N_TCLK;clk_id++){
	for(int str_id=0;str_id<N_STRP;str_id++){
//****	      if(tpc_data[0][clk_id][str_id]>0) h_anode->SetBinContent(str_id,clk_id,tpc_data[0][clk_id][str_id]);
//****	      if(tpc_data[1][clk_id][str_id]>0) h_cathode->SetBinContent(str_id,clk_id,tpc_data[1][clk_id][str_id]);
	  if(tpc_data[0][clk_id][str_id]>0) h_anode->SetBinContent(str_id,clk_id,1);
	  if(tpc_data[1][clk_id][str_id]>0) h_cathode->SetBinContent(str_id,clk_id,1);
	  h_raw_anode->SetBinContent(str_id,clk_id,raw_wave[0][clk_id][str_id]);
	  h_raw_cathode->SetBinContent(str_id,clk_id,raw_wave[1][clk_id][str_id]);
	}
      }
#ifdef DRAW_CANVAS
      c1->cd();
      c1->Clear();
      h_raw_anode->Draw("colz");
#ifdef DRAW_ARROW
      for(Int_t i_alpha=0;i_alpha<nAlpha;i_alpha++){
	ar_anode[i_alpha].Draw();
      }
#endif
#endif
#ifdef DRAW_CANVAS
      c1->Update();
      c2->cd();
      c2->Clear();
      h_raw_cathode->Draw("colz");
#ifdef DRAW_ARROW
      for(Int_t i_alpha=0;i_alpha<nAlpha;i_alpha++){
	ar_cathode[i_alpha].Draw();
      }
#endif
#endif
#ifdef DRAW_CANVAS
      c2->Update();
      c3->cd();
      c3->Clear();
      h_anode->Draw("box");
#ifdef DRAW_ARROW
      for(Int_t i_alpha=0;i_alpha<nAlpha;i_alpha++){
	ar_anode[i_alpha].Draw();
      }
#endif
#endif
#ifdef DRAW_CANVAS
      c3->Update();
      c4->cd();
      c4->Clear();
      h_cathode->Draw("box");
#ifdef DRAW_ARROW
      for(Int_t i_alpha=0;i_alpha<nAlpha;i_alpha++){
	ar_cathode[i_alpha].Draw();
      }
#endif
#endif
#ifdef DRAW_CANVAS
      c4->Update();
      std::cout << "update" << std::endl;
      h_raw_anode->Write();
      h_raw_cathode->Write();
      h_anode->Write();
      h_cathode->Write();
      if(getchar() == 'q') return 0;
#endif
#ifdef SAVE_FILE
//****      tree->Fill();
      for(int a=0;a<N_AC;a++){
	for(int b=0;b<N_TCLK;b++){
	  for(int c=0;c<N_STRP;c++){
	    fprintf(file_tpc,"%d ",tpc_data[a][b][c]);
	  }
	}
      }
      fprintf(file_tpc,"\n");
      for(int i_alpha=0;i_alpha<nAlpha;i_alpha++){
	for(int a=0;a<2;a++){
	  fprintf(file_para,"%f ",recoil_track_a_x[i_alpha][a]);
	  fprintf(file_para,"%f ",recoil_track_a_y[i_alpha][a]);
	  fprintf(file_para,"%f ",recoil_track_c_x[i_alpha][a]);
	  fprintf(file_para,"%f ",recoil_track_c_y[i_alpha][a]);
	}
	fprintf(file_para,"%f ",e3[i_alpha]);
	fprintf(file_para,"%f ",theta3_deg[i_alpha]);
	fprintf(file_para,"%f ",phi3_deg[i_alpha]);
	fprintf(file_para,"%f ",dr[i_alpha]);
      }
      fprintf(file_para,"\n");
#endif
      h_anode->Reset();
      h_cathode->Reset();
      h_raw_anode->Reset();
      h_raw_cathode->Reset();
      
      
      ii++;
    }  // end of if(alpha_inside...)
  }  // end of ii loop
#ifdef SAVE_FILE
//****  tree->Write();
//****  tfile->Close();
//****  delete tfile;
  fclose(file_tpc);
  fclose(file_para);
#endif
  
  //      printf("theta=%.2f, e3=%.2f, eff=%.4f\n", theta3_deg, e3, stop_eff);
  stop_eff = stop_cnt/((double)total_cnt);
  
  //**  }  // end of theta3 loop
  
  /* save dat file */
//****  for(int ii=0;ii<filesize;ii++){
//****    fclose(fp_t[ii]);
//****    fclose(fp_r[ii]);
//****  }
#ifndef SAVE_FILE
  app.Run();
#endif
  
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

