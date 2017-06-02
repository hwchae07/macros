#ifndef ANAHOD_H
#define ANAHOD_H

#include "AnaModule.H"

#include "TArtSAMURAIParameters.hh"
#include "TArtCalibHODPla.hh"

class AnaHOD : public AnaModule {

public:

  AnaHOD();
  virtual ~AnaHOD();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtSAMURAIParameters* fHODParameters;
  
  TArtCalibHODPla* fCalibHODPla;

  // for trees
  Int_t hodNum; // num Pla
  Int_t* hodID;//[24]; // ID
  Double_t* hodPlaQU;//[24]; // Cal
  Double_t* hodPlaQD;//[24]; // Cal
  Double_t* hodPlaTU;//[24]; // Slw
  Double_t* hodPlaTD;//[24]; // Slw
  Double_t* hodPlaTA;//[24]; // Slw
  Double_t *hodTURaw;
  Double_t *hodTDRaw;

  Double_t *hodQURaw;
  Double_t *hodQDRaw;
  
  Double_t *hodQ;
  
  // Tof offset
  Double_t parTOF[25]; // [0]: total offset, 1 ~ 24: each
};

#endif
