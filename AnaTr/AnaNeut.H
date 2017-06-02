#ifndef ANANEUT_H
#define ANANEUT_H

#include "AnaModule.H"

#include "TArtSAMURAIParameters.hh"
#include "TArtCalibNEBULA.hh"
#include "TArtCalibNeuLAND.hh"
#include "TArtCalibNeuLANDVETO.hh"

class AnaNeut : public AnaModule {

public:

  AnaNeut();
  virtual ~AnaNeut();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtSAMURAIParameters* fNeutParameters;
  
  TArtCalibNEBULA* fCalibNEBULA;
  TArtCalibNeuLAND* fCalibNeuLAND;
  TArtCalibNeuLANDVETO* fCalibNeuLANDVETO;

  // for parameters
  Double_t* tofOffsetNEB;
  Double_t* tofOffsetNEU;

  // for trees
  Int_t neutNum;
  Int_t neulandNum;
  Int_t nebulaNum;
  Int_t neulandMult[9];
  Int_t nebulaMult[6];

  Int_t* neutID;//[144];
  Double_t* neutQU;//[144]; // Cal
  Double_t* neutQD;//[144]; // Cal
  Double_t* neutTU;//[144]; // Slw
  Double_t* neutTD;//[144]; // Slw
  Double_t* neutTA;//[144]; // Slw
  Double_t* neutX;//[144];
  Double_t* neutY;//[144];
  Double_t* neutZ;//[144];  
  Bool_t* neutNEB;
  Bool_t* neutVETO;

};

#endif
