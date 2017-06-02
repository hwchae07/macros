#ifndef ANABDC_H
#define ANABDC_H

#include "AnaModule.H"

#include "TArtSAMURAIParameters.hh"
#include "TArtCalibBDC1Hit.hh"
#include "TArtCalibBDC1Track.hh"
#include "TArtCalibBDC2Hit.hh"
#include "TArtCalibBDC2Track.hh"

class AnaBDC : public AnaModule {

public:

  AnaBDC();
  virtual ~AnaBDC();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  void LoadTDCDistribution(const char* filename);
  TArtSAMURAIParameters* fBDCParameters;
  
  TArtCalibBDC1Hit*   fCalibBDC1Hit;
  TArtCalibBDC1Track* fCalibBDC1Track;
  TArtCalibBDC2Hit*   fCalibBDC2Hit;
  TArtCalibBDC2Track* fCalibBDC2Track;


  // for trees
  Double_t bdc1X, bdc1Y, bdc1A, bdc1B, bdc1CX, bdc1CY;
  Double_t bdc2X, bdc2Y, bdc2A, bdc2B, bdc2CX, bdc2CY;
};

#endif
