#include "12C_3alpha.hpp"
#include <random>

void decay12C(TLorentzVector decay_alpha[], int particle[])
{
  TRandom *rndm = new TRandom();
  std::random_device seed_gen;

  Double_t E_beam = T_beam+M_n;
  Double_t Pz_beam = TMath::Sqrt(E_beam*E_beam-M_n*M_n);
  TLorentzVector beam(0.0,0.0,Pz_beam,E_beam);
  TLorentzVector target(0.0,0.0,0.0,M_10C);
  TLorentzVector W = beam+target;
  Double_t masses1[2] = {M_alpha,M_10C};
  particle[0] = 0; // 0:4He 1:10C 2:12C
  particle[1] = 1;
//  Doutle_t masses1[2] = {M_12C, M_10C};
//  particle[0] = 2;
//  particle[1] = 1;

  TGenPhaseSpace *event = new TGenPhaseSpace();

  Double_t weight;
  Double_t uniform_rndm;

  rndm->SetSeed(seed_gen());
  event->SetDecay(W,2,masses1);
  do{
    weight = event->Generate();
    uniform_rndm = rndm->Uniform(0.,event->GetWtMax());
  }while(uniform_rndm > weight);
  for(Int_t n=0;n<2;n++){
    decay_alpha[n] = *event->GetDecay(n);
  }

  delete rndm;
  delete event;
}

int decay12C(std::vector<TLorentzVector> *decay_alpha,int flag)
{
  TRandom *rndm = new TRandom();
  std::random_device seed_gen;

  Double_t E_beam = T_beam+M_n; //total energy of incident neutron [GeV]
  Double_t Pz_beam = TMath::Sqrt(E_beam*E_beam-M_n*M_n); //Pz of neutron [GeV]
  TLorentzVector beam(0.0,0.0,Pz_beam,E_beam);
  TLorentzVector target(0.0,0.0,0.0,M_12C);
  TLorentzVector W = beam+target;
  Double_t masses1[2] = {M_12C,M_n};

  TGenPhaseSpace *event = new TGenPhaseSpace();
  switch(flag){
  case 0:
    {
      masses1[0] = masses1[0] + 7.644; // 0^+_2
    break;
    }
  case 1:
    {
      masses1[0] = masses1[0] + 9.870; // 2^+_2
    break;
    }
  default:
    {
      masses1[0] = masses1[0] + Ex_12C; // default value
    break;
    }
  }
  Double_t masses2[nAlpha] = {M_alpha,M_alpha};
  Double_t E_alpha[3];

  TLorentzVector *sca1[2];
  Double_t weight;
  Double_t uniform_rndm;

  rndm->SetSeed(seed_gen());

  event->SetDecay(W,2,masses1);
  do{
    weight = event->Generate();
    uniform_rndm = rndm->Uniform(0.,event->GetWtMax());
  }while(uniform_rndm > weight);
  for(Int_t n=0;n<2;n++){
    sca1[n] = event->GetDecay(n);
  }
  
  event->SetDecay(*sca1[0],nAlpha,masses2);
  do{
    weight = event->Generate();
    uniform_rndm = rndm->Uniform(0.,event->GetWtMax());
  }while(uniform_rndm > weight);

  for(Int_t n=0;n<nAlpha;n++){
    decay_alpha->push_back(*event->GetDecay(n));
  }

  // clear pointer
  delete rndm;
  delete event;
  
  return nAlpha;
}

int BG(std::vector<TLorentzVector> *back_ground,int flag)
{
  TRandom *rndm = new TRandom();
  std::random_device seed_gen;

  Double_t E_beam = T_beam+M_n;
  Double_t Pz_beam = TMath::Sqrt(E_beam*E_beam-M_n*M_n);
  TLorentzVector beam(0.0,0.0,Pz_beam,E_beam);
  Int_t n_particle;
  Double_t masses[2];
  TLorentzVector target;

  switch(flag){
  case 0:
    {
      masses[0] = M_12C;
      masses[1] = M_n;
      target = TLorentzVector(0.0,0.0,0.0,M_12C);
      n_particle = 2;
      break;
    }
  case 1:
    {
      masses[0] = M_p;
      masses[1] = M_n;
      target = TLorentzVector(0.0,0.0,0.0,M_p);
      n_particle = 2;
      break;
    }
  default:
    {
      return 0;
    }
  }

  TLorentzVector W = beam+target;

  TGenPhaseSpace *event = new TGenPhaseSpace();
  Double_t weight;
  Double_t uniform_rndm;

  rndm->SetSeed(seed_gen());

  event->SetDecay(W,n_particle,masses);
  do{
    weight = event->Generate();
    uniform_rndm = rndm->Uniform(0.,event->GetWtMax());
  }while(uniform_rndm > weight);
  for(Int_t n=0;n<n_particle;n++){
    back_ground->push_back(*event->GetDecay(n));
  }

  delete rndm;
  delete event;

  return n_particle;
}
