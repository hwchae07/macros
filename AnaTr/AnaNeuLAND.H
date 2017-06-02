#ifndef ANANEULAND_H
#define ANANEULAND_H

#include "AnaModule.H"

#include "TArtSAMURAIParameters.hh"
#include "TArtCalibNeuLAND.hh"
#include "TArtCalibNeuLANDVETO.hh"

class AnaNeuLAND : public AnaModule {

public:

  AnaNeuLAND();
  virtual ~AnaNeuLAND();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtSAMURAIParameters* fNeuLANDParameters;
  
  TArtCalibNeuLAND* fCalibNeuLAND;
  TArtCalibNeuLANDVETO* fCalibNeuLANDVETO;

  // for trees
  Int_t neulandNum;
  Int_t neulandMult[9];

  //  Double_t neulandTRef; // tref
  Int_t* neulandID;//[144];
  Double_t* neulandQU;//[144]; // Cal
  Double_t* neulandQD;//[144]; // Cal
  Double_t* neulandTU;//[144]; // Slw?
  Double_t* neulandTD;//[144]; // Slw?
  Double_t* neulandTA;//[144]; // Slw?
  Double_t* neulandX;//[144];
  Double_t* neulandY;//[144];
  Double_t* neulandZ;//[144];

  Double_t* neulandTURaw;//[144];
  Double_t* neulandTDRaw;//[144];

  //Double_t* tofOffset;
};

#endif
