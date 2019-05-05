#ifndef _12_3ALPHA_HPP_
#define _12_3ALPHA_HPP_

#include <TROOT.h>
#include <TLorentzVector.h>
#include <TGenPhaseSpace.h>
#include <TRandom.h>
#include <vector>
#include "para.h"

const Double_t GeV = 1.;
const Double_t MeV = 1000.;
const Double_t keV = 1000000.;
const Double_t T_beam = beam_ene/MeV; //kinetic energy of incident neutron [GeV]
const Double_t M_alpha = Mass_4He/MeV; //alpha mass [GeV]
const Double_t M_12C = Mass_12C/MeV; //12C mass [GeV]
const Double_t Ex_12C = 7.654/MeV; //12C excited energy [GeV] (0^+_2)
const Double_t M_n = Mass_n/MeV; //n mass [GeV]
const Double_t M_p = Mass_p/MeV; //p mass [GeV]

void decay12C(TLorentzVector decay_alpha[]);
int decay12C(std::vector<TLorentzVector> *decay_alpha,int flag);
int BG(std::vector<TLorentzVector> *back_ground,int flag);

#endif
