#ifndef ANADALI_H
#define ANADALI_H

#include "AnaModule.H"

#include "TArtDALIParameters.hh"
#include "TArtCalibDALI.hh"

class AnaDALI : public AnaModule {

public:

  AnaDALI();
  virtual ~AnaDALI();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtDALIParameters* fDALIParameters;
  
  TArtCalibDALI* fCalibDALI;

  // for trees
  Int_t daliNum;

  Int_t* daliID;//[144];
  Double_t* daliTheta;//[144]; // Cal
  Double_t* daliE;//[144]; // Cal
  Double_t* daliT;//[144]; // Cal
};

#endif
