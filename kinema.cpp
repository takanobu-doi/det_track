#include "kinema.h"

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
#include <TCanvas.h>
#include <TStyle.h>

Double_t calc_Ex(Double_t m1, Double_t m2, Double_t m3, Double_t m4,
                 Double_t K1, Double_t K3, Double_t theta3){

  Double_t E1, E3;
  Double_t p1, p3;
  Double_t s,t,u;
  Double_t total_m4;
  Double_t ex4;

  E1 = m1+K1;
  E3 = m3+K3;

  p1 = sqrt(E1*E1-m1*m1);
  p3 = sqrt(E3*E3-m3*m3);

  s=m1*m1+m2*m2+2*m2*E1;
  t=m1*m1+m3*m3+2*(p1*p3*cos(theta3)-E1*E3);
  u=m2*m2+m3*m3-2*m2*E3;

  total_m4=sqrt(s+t+u-m1*m1-m2*m2-m3*m3);
  ex4=total_m4-m4;

  return ex4;
}

Double_t calccmang(Double_t m1, Double_t m2, Double_t m3, Double_t m4,
		   Double_t beamene, Double_t ex, Double_t labang_r, Double_t e3){
  
  double p1,p2; // momentum in LAB
  double p1c,p2c,p3c,p4c; //momentum in CM
  double E1,E2; //total energy in LAB
  double E1c,E2c,E3c,E4c; //total energy in CM

  double Q_val; // Q-value of the reaction
  double W; // total CM energy
  double gamma_cm;
  double beta_cm;
  double g2c; //gamma factor in CM
  double b1c,b2c,b3c,b4c; //velocity in CM
  double delta23,delta24; //velocity ratio
  
  E1=m1+beamene;
  E2=m2;
  m4+=ex;
  
  Q_val = m3+m4-m1-m2;
  p1=sqrt(E1*E1-m1*m1);
  p2=0;
  W=sqrt(m1*m1+m2*m2+2*m2*E1);
  gamma_cm=(E1+m2)/W;
  beta_cm=p1/(E1+m2);
  p1c=1.0/(2*W)*sqrt((W*W-(m1+m2)*(m1+m2))*(W*W-(m1-m2)*(m1-m2)));
  p2c=p1c;
  p3c=1.0/(2*W)*sqrt((W*W-(m3+m4)*(m3+m4))*(W*W-(m3-m4)*(m3-m4)));
  p4c=p3c;
  E1c=sqrt(m1*m1+p1c*p1c);
  E2c=sqrt(m2*m2+p2c*p2c);
  E3c=sqrt(m3*m3+p3c*p3c);
  E4c=sqrt(m4*m4+p4c*p4c);
  g2c=E2c/m2;
  b1c=p1c/E1c;
  b2c=p2c/E2c;
  b3c=p3c/E3c;
  b4c=p4c/E4c;
  delta23=b2c/b3c;
  delta24=b2c/b4c;

  double firstterm,secondterm;
  int n=0;
  double xsq,tmpang;
  double cmcos[2];

  ///////////////////////////////
  double gamma2 = g2c;
  double delta = delta23;

  xsq=pow(gamma2*tan(labang_r),2.0);

  firstterm=-delta*xsq/(1.0+xsq);
  secondterm=sqrt((1.0-delta*delta)*xsq+1)/(1.0+xsq);

  if( ((1.0-delta*delta)*xsq+1) < 0) printf("CM ang cal error!\n");
  
  cmcos[0]=firstterm+secondterm;
  cmcos[1]=firstterm-secondterm;

  double ang[2], ene[2];
  double result_ang;
  ang[0]=-100;
  ang[1]=-100;  
  ene[0]=-100;
  ene[1]=-100;

  
  if(fabs(cmcos[0])< 1 || cmcos[0]==1) {
    tmpang=acos(cmcos[0]);
    if((cos(tmpang)+delta)*cos(labang_r)>0){
      ang[0]=tmpang;
      ene[0]=gamma2*(b2c*p3c*cos(tmpang)+E3c)-m3;
      n++;
    }
  }

  if(((fabs(cmcos[1])<=1) && (fabs(secondterm)>10e-10)) || (cmcos[1]==-1)) {
    tmpang=acos(cmcos[1]);
    if((cos(tmpang)+delta)*cos(labang_r)>0){
      ang[1]=tmpang;
      ene[1]=gamma2*(b2c*p3c*cos(tmpang)+E3c)-m3;
      n++;
    }
  }

  if(n==1){
    if(ang[0]>-100) result_ang=180-ang[0]*TMath::RadToDeg();
    if(ang[1]>-100) result_ang=180-ang[1]*TMath::RadToDeg();    
  }

  if(n==2){
    if(fabs(e3-ene[0]) <  fabs(e3-ene[1])) result_ang=180-ang[0]*TMath::RadToDeg();
    if(fabs(e3-ene[0]) >= fabs(e3-ene[1])) result_ang=180-ang[1]*TMath::RadToDeg();    
  }
  
  return(result_ang);

}
