#ifndef ANAFDC_H
#define ANAFDC_H

#include "AnaModule.H"

#include "TArtSAMURAIParameters.hh"
#include "TArtCalibFDC1Hit.hh"
#include "TArtCalibFDC1Track.hh"
#include "TArtCalibFDC2Hit.hh"
#include "TArtCalibFDC2Track.hh"

class AnaFDC : public AnaModule {

public:

  AnaFDC();
  virtual ~AnaFDC();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  void LoadTDCDistribution(const char* filename);
  TArtSAMURAIParameters* fFDCParameters;
  TArtCalibFDC1Hit*   fCalibFDC1Hit;
  TArtCalibFDC1Track* fCalibFDC1Track;
  TArtCalibFDC2Hit*   fCalibFDC2Hit;
  TArtCalibFDC2Track* fCalibFDC2Track;


  // for trees
  Double_t fdc1X, fdc1Y, fdc1A, fdc1B, fdc1C;
  Double_t fdc2X, fdc2Y, fdc2A, fdc2B, fdc2C;

  Double_t fdc1ex,fdc1ey,fdc1ez;
};
#endif
