#include <iostream>
#include <cstdlib>
#include <string>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include "database.hpp"
#include "kinema_doi.hpp"

kinema::kinema()
{
  Ex2 = 0;
  thr = 0;
  beta_c = 0;
}

kinema::kinema(int argc,char *argv[])
{
  Ex2 = 0;
  thr = 0;
  beta_c = 0;
  setenv(argc,argv);
}

int kinema::calc()
{
  double a,b,c;
  double M;
  double E0;

  setparameter();

  M = sqrt(pow(p[0].get_mass(),2)+pow(p[1].get_mass(),2)+2*p[1].get_E()*p[0].get_E());
  E0 = (M*M+pow(p[2].get_mass(),2)-pow(p[3].get_mass(),2))/(2*M);

  a = 1-pow(beta_c*cos(M_PI*thr/180.),2);
  b = -2*E0*sqrt(1-beta_c*beta_c);
  c = E0*E0*(1-beta_c*beta_c)+pow(beta_c*p[2].get_mass()*cos(M_PI*thr/180.),2);

  if(thr<90.){
    p[2].set_E((-b+sqrt(fabs(b*b-4*a*c)))/(2*a)-p[2].get_mass());
  }else{
    p[2].set_E((-b-sqrt(fabs(b*b-4*a*c)))/(2*a)-p[2].get_mass());
  }
  P = sqrt(pow(p[2].get_E(),2)-pow(p[2].get_mass(),2));
  p0 = sqrt((pow((M*M+pow(p[2].get_mass(),2)-pow(p[3].get_mass(),2))/(2*M),2)-pow(p[2].get_mass(),2)));

  //  p[3].set_E(sqrt(pow(p[3].get_mass(),2)+pow((p0*cos(M_PI*getthr3_C()/180.)+(M*M+pow(p[3].get_mass(),2)-pow(p[2].get_mass(),2))/(2*M)*beta_c)/sqrt(1-beta_c*beta_c),2)+pow(p0*sin(M_PI*getthr3_C()/180.),2))-p[3].get_mass());
  p[3].set_E(sqrt(pow(p[3].get_mass(),2)+pow(p0,2))/sqrt(1-pow(beta_c,2))-p0*cos(getthr3_C()*M_PI/180.)*beta_c/sqrt(1-pow(beta_c,2))-p[3].get_mass());

  
  if(getthr3()>getthr3_max()){
    return 1;
  }else{
    return 0;
  }
}

int kinema::calc(double E_inc,double thr_L,double ex3,double ex4)
{
  double a,b,c;
  double M;
  double E0;

  setparameter(E_inc,thr_L,ex3,ex4);

  M = sqrt(pow(p[0].get_mass(),2)+pow(p[1].get_mass(),2)+2*p[1].get_E()*p[0].get_E());
  E0 = (M*M+pow(p[2].get_mass(),2)-pow(p[3].get_mass(),2))/(2*M);

  a = 1-pow(beta_c*cos(M_PI*thr/180.),2);
  b = -2*E0*sqrt(1-beta_c*beta_c);
  c = E0*E0*(1-beta_c*beta_c)+pow(beta_c*p[2].get_mass()*cos(M_PI*thr/180.),2);

  if(thr<90.){
    p[2].set_E((-b+sqrt(b*b-4*a*c))/(2*a)-p[2].get_mass());
  }else{
    p[2].set_E((-b-sqrt(b*b-4*a*c))/(2*a)-p[2].get_mass());
  }
  P = sqrt(pow(p[2].get_E(),2)-pow(p[2].get_mass(),2));
  p0 = sqrt((pow((M*M+pow(p[2].get_mass(),2)-pow(p[3].get_mass(),2))/(2*M),2)-pow(p[2].get_mass(),2)));

  p[3].set_E(sqrt(pow(p[3].get_mass(),2)+pow(p0,2))/sqrt(1-pow(beta_c,2))-p0*cos(getthr3_C()*M_PI/180.)*beta_c/sqrt(1-pow(beta_c,2))-p[3].get_mass());
  
  if(getthr3()>getthr3_max()){
    return 1;
  }else{
    return 0;
  }
}

int kinema::ex4toei()
{
  double X,Y;
  double E0;
  double E2,Ex;
  int i;
  char s1[100];

  fgets(s1,100,stdin);
  i = sscanf(s1,"%lf %lf %lf",&E2,&thr,&Ex);
  if(i==1){
    Ex = 0;
    thr = 0;
    i = 3;
  }if(i==2){
    Ex = 0;
    i = 3;
  }if(i!=3){
    cout << "Illigal format." << endl;
    exit(EXIT_SUCCESS);
  }
  p[1].set_E(0);
  p[2].set_E(E2);
  p[3].set_Ex(Ex);

  X = (-2*p[1].get_E()*p[2].get_E()+pow(p[1].get_E(),2)-pow(p[3].get_mass(),2)+pow(p[0].get_mass(),2)+pow(p[2].get_mass(),2))/(2*(p[1].get_E()-p[2].get_E()));
  Y = (pow(p[2].get_E(),2)-pow(p[2].get_mass(),2))*pow(cos(M_PI*thr/180.),2)/pow(p[1].get_E()-p[2].get_E(),2);

  E0 = (-X-sqrt(X*X-(1-Y)*(X*X+pow(p[0].get_mass(),2)*Y)))/(1-Y);
  
  if(E0>p[0].get_mass()){
     setE_inc(E0-p[0].get_mass());
    return 0;
  }else{
    return 1;
  }
}

double kinema::get_gfactor()
{
  double u_cm,gamma,E0;

  gamma = pow(1-beta_c*beta_c,-0.5);
  E0 = gamma*p[2].get_E()-beta_c*gamma*sqrt(pow(p[2].get_E(),2)-pow(p[2].get_mass(),2))*cos(thr*M_PI/180.);
  u_cm = sqrt(1-pow(p[2].get_mass(),2)/(E0*E0));

  return gamma*(1+cos(getthr3_C()*M_PI/180.)*beta_c/u_cm)/pow(pow(sin(getthr3_C()*M_PI/180.),2)+gamma*gamma*pow(cos(getthr3_C()*M_PI/180.)+beta_c/u_cm,2),3./2.);
}

void kinema::setenv(int argc,char *argv[])
{
  if(argc==6){
    setparticles(argv[2],argv[3],argv[4],argv[5]);
    setE_inc(strtod(argv[1],NULL));
  }else{
    cout << "Illigal format." << endl;
    exit(EXIT_FAILURE);
  }
  return;
}
/*
void kinema::setparticles()
{
  cout << endl;
  cout << "choose particle (0:project 1:target 2:scatter 3:recoil)";
  cout << "ex:12C, p, etc..." << endl;
  for(int i=0;i<4;i++){
    cout << i+1;
    while(1){
      if(set.read_stdin()==0){
	m[i] = set.get_A()*amu + set.get_delta(set.get_A(),set.get_name());
	a[i] = set.get_A();
	n[i] = set.get_name();
	break;
      }
      cout << "please again" << endl << i+1;
    }
  }
}
*/
void kinema::setparticles(char *n0,char *n1,char *n2,char *n3)
{
  char *ss[4]={n0,n1,n2,n3};
  for(int i=0;i<4;i++){
    if(set.read_comd(ss[i])==0){
      p[i].set_delta(set.get_delta(set.get_A(),set.get_name()));
      p[i].set_N(set.get_N(set.get_A(),set.get_name()));
      p[i].set_Z(set.get_Z(set.get_A(),set.get_name()));
      p[i].set_name(set.get_name());
    }else{
      cerr << "I cannot find the particle that you wrote." << endl;
      exit(EXIT_FAILURE);
    }
  }
  if(p[0].get_N()+p[1].get_N()-p[2].get_N()-p[3].get_N()!=0){
    if(p[0].get_Z()+p[1].get_Z()-p[2].get_Z()-p[3].get_Z()!=0){
      cerr << "This reaction cannot realize." << endl;
      exit(EXIT_FAILURE);
    }
  }
}  

void kinema::showparticles()
{
  for(int i=0;i<4;i++){
    cout << "#" << i+1 << ": " << setw(3) << right << p[i].get_A();
    cout << setw(2) << left << p[i].get_name() << " m=";
    cout << setw(10) << right << fixed << setprecision(3) << p[i].get_mass();
    cout << "(MeV)" << " delta = " << p[i].get_delta() << "(MeV)" << endl;
  }
}

void kinema::setparameter(double E_inc,double thr_L,double ex3,double ex4)
{
  Ex1 = ex3;
  Ex2 = ex4;
  thr = thr_L;
  p[0].set_E(E_inc);
  p[1].set_E(0);
  p[2].set_Ex(Ex1);
  p[3].set_Ex(Ex2);
  beta_c = sqrt(p[0].get_E()*p[0].get_E()-p[0].get_mass()*p[0].get_mass())/(p[0].get_E()+p[1].get_E());
}

void kinema::setparameter(double E)
{
  double thr,ex3,ex4;
  char s1[100];
  int i;
  //  cout << "thr_L Ex4" << endl;
  //  cin >> thr >> Ex;
  fgets(s1,100,stdin);
  i = sscanf(s1,"%lf %lf %lf",&thr,&ex4,&ex3);
  if(i==1){
    ex3 = 0;
    ex4 = 0;
    i = 3;
  }
  else if(i==2){
    ex3 = 0;
    i = 3;
  }
  if(i!=3){
    cout << "Illigal format." << endl;
    exit(EXIT_SUCCESS);
  }
  setparameter(E,thr,ex3,ex4);
}

void kinema::setparameter()
{
  double K;
  if(p[0].get_K()==0){
    cout << "E_inc" << endl;
    cin >> K;
    p[0].set_E(K);
  }
  setparameter(p[0].get_K());
}

void kinema::setE_inc()
{
  double K;
  cout << "E_inc" << endl;
  cin >> K;
  p[0].set_E(K);
}

void kinema::setE_inc(double E_inc)
{
  p[0].set_E(E_inc);
}

void kinema::showparameter()
{
  if(p[0].get_E()!=0){
    cout << "# E_inc = " << p[0].get_K() << endl;
  }
}


double kinema::getE_L(int i)
{
  return p[i].get_K();
}

double kinema::getthr3()
{
  return thr;
}

double kinema::getthr3_C()
{
  double M,E0,thr_c,gamma;
  gamma = pow(1-beta_c*beta_c,-0.5);
  M = sqrt(pow(p[0].get_mass(),2)+pow(p[1].get_mass(),2)+2*p[1].get_E()*p[0].get_E());
  E0 = (M*M+pow(p[2].get_mass(),2)-pow(p[3].get_mass(),2))/(2*M);
  /*
  if(thr==0){
    return 0.;
    }
  thr_c = atan((sqrt(1-beta_c*beta_c)*p*cos(M_PI*thr/180.)-(E0*beta_c))/p*sin(M_PI*thr/180.))*180/M_PI;
  */
  thr_c = atan(sqrt(pow(p[2].get_E(),2)-pow(p[2].get_mass(),2))*sin(thr*M_PI/180.)/(-beta_c*gamma*(p[2].get_E())+gamma*sqrt(pow(p[2].get_E(),2)-pow(p[2].get_mass(),2))*cos(thr*M_PI/180.)))*180./M_PI;
  if(thr_c>=0){
    return thr_c;
  }else{
    return thr_c+180.;
  }
}

double kinema::getthr3_max()
{
  double theta,M;

  M = sqrt(pow(p[0].get_mass(),2)+pow(p[1].get_mass(),2)+2*p[1].get_E()*p[0].get_E());
  if(p[0].get_mass()==p[1].get_mass()&&
     p[1].get_mass()==p[2].get_mass()&&
     p[2].get_mass()==p[3].get_mass()){
    return 90.;
  }
  if(beta_c*(M*M+pow(p[2].get_mass(),2)-pow(p[3].get_mass(),2))/(2*M)<p0){
    return 180.;
  }
  theta = asin(p0*sqrt(1-beta_c*beta_c)/(p[2].get_mass()*beta_c))*180./M_PI;
  if(theta>=0){
    return theta;
  }else{
    return -theta;
  }
}

double kinema::getthr4()
{
  double momentum[4] = {0,0,0,0};
  double cos_psi,sin_psi,tan_psi,cot_psi;
  double psi;
  /*
  for(int i=0;i<4;i++){
    momentum[i] = sqrt(pow(p[i].get_E(),2)-pow(p[i].get_mass(),2));
  }
  sin_psi = momentum[2]*sin(thr*M_PI/180.)/momentum[3];
  psi = asin(sin_psi)*180./M_PI;
  */
  double p1,p3;
  p1 = sqrt(pow(p[0].get_E(),2)-pow(p[0].get_mass(),2));
  p3 = sqrt(pow(p[2].get_E(),2)-pow(p[2].get_mass(),2));
  cot_psi = (p1-p3*cos(thr*M_PI/180.))/(p3*sin(thr*M_PI/180.));
  psi = 90-atan(cot_psi)*180./M_PI;
  
  return psi;
}

double kinema::getEx1()
{
  return Ex1;
}

double kinema::getEx2()
{
  return Ex2;
}

double kinema::getmass(int i)
{
  return p[i].get_mass();
}

double kinema::getp()
{
  return P;
}

double kinema::Qvalue()
{
  return p[0].get_delta()+p[1].get_delta()-p[2].get_delta()-p[3].get_delta();
}
