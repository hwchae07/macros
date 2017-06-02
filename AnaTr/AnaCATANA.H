#ifndef ANACATANA_H
#define ANACATANA_H

#include "AnaModule.H"

#include "TArtCATANAParameters.hh"
#include "TArtCalibCATANA.hh"

class AnaCATANA : public AnaModule {

public:

  AnaCATANA();
  virtual ~AnaCATANA();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtCATANAParameters* fCATANAParameters;
  
  TArtCalibCATANA* fCalibCATANA;

  // for trees
  Int_t catanaNum;

  Int_t* catanaID;//[144];
  Double_t* catanaTheta;//[144]; // Cal
  Double_t* catanaE;//[144]; // Cal
  Double_t* catanaQ;//[144]; // Cal
  Double_t* catanaT;//[144]; // Cal
  Double_t *catanaRawTDC;
  Double_t *catanaRawTRef;
};

#endif
