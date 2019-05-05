#ifndef _KINEMA_HPP_
#define _KINEMA_HPP_

#include "database.hpp"
#include "nuclear.hpp"

class kinema
{
private:
  database set;
  /*
  double m[4];
  string n[4];
  int a[4];
  */
  nucl p[4];
  double Ex1;
  double Ex2;
  double thr;
  double beta_c;
  double p0,P;
public:
  kinema();
  kinema(int argc,char *argv[]);
  void setenv(int argc,char *argv[]);
  //  void setparticles();
  void setparticles(char*,char*,char*,char*);
  void setE_inc();
  void setE_inc(double);
  void showparticles();
  void setparameter();//(if faild return 1)
  void setparameter(double);//E_inc (if faild return 1)
  void setparameter(double,double,double,double);//E_inc,thr_L,Ex
  void showparameter();
  int calc();
  int calc(double E_inc,double thr_L,double ex3,double ex4);
  int ex4toei();
  double get_gfactor();
  double getE_L(int);
  double getthr3();
  double getthr3_C();//計算にバグが残っている(丸め誤差?)
  double getthr3_max();
  double getthr4();
  double getEx1();
  double getEx2();
  double getmass(int);
  double getp();
  double Qvalue();
};

#endif
