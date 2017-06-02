#ifndef ANAPPAC_H
#define ANAPPAC_H

#include "AnaModule.H"

#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPPAC.hh"

class AnaPPAC : public AnaModule {

public:

  AnaPPAC();
  virtual ~AnaPPAC();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtBigRIPSParameters* fPPACParameters;
  
  TArtCalibPPAC* fCalibPPAC;


  // for trees
  Double_t f5PPAC1AX, f5PPAC1AY;
  Double_t f5PPAC1BX, f5PPAC1BY;
  Double_t f5PPAC2AX, f5PPAC2AY;
  Double_t f5PPAC2BX, f5PPAC2BY;

  Double_t f5PPAC1AQA,f5PPAC1ATA;
  Double_t f5PPAC1BQA,f5PPAC1BTA;
  Double_t f5PPAC2AQA,f5PPAC2ATA;
  Double_t f5PPAC2BQA,f5PPAC2BTA;
};

#endif
