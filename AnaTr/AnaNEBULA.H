#ifndef ANANEBULA_H
#define ANANEBULA_H

#include "AnaModule.H"

#include "TArtSAMURAIParameters.hh"
#include "TArtCalibNEBULA.hh"

class AnaNEBULA : public AnaModule {

public:

  AnaNEBULA();
  virtual ~AnaNEBULA();

  void InitParameter();
  void InitDetector();
  void Analysis();

  void SetTree();

  void DeleteAll();

private:

  TArtSAMURAIParameters* fNEBULAParameters;

  TArtCalibNEBULA* fCalibNEBULA;

  // for parameters
  //Double_t* tofOffset;
  double **parTU;
  double **parTD;

  double* par_ycal0;
  double* par_ycal1;
  double* TAveOff;

  double* tdcAlignU;
  double* tdcAlignD;

  // for trees
  Int_t nebulaNum;
  Int_t nebulaMult[6];

  Int_t* nebulaID;//[144];
  Double_t* nebulaQU;//[144]; // Cal
  Double_t* nebulaQD;//[144]; // Cal
  Double_t* nebulaTU;//[144]; // Slw
  Double_t* nebulaTD;//[144]; // Slw
  Double_t* nebulaTA;//[144]; // Slw
  Double_t* nebulaX;//[144];
  Double_t* nebulaY;//[144];
  //Double_t* nebulaY2;
  Double_t* nebulaZ;//[144];
  Double_t* nebulaTURaw;
  Double_t* nebulaTDRaw;

  Double_t* nebulaTU2;
  Double_t* nebulaTD2;
  Double_t* nebulaTA2;
  Double_t* nebulaDT;

  Double_t* nebulaQURaw;
  Double_t* nebulaQDRaw;
  Double_t* nebulaQUPed;
  Double_t* nebulaQDPed;
  Double_t* nebulaQPed;

  Double_t* nebulaQA;
};

#endif
