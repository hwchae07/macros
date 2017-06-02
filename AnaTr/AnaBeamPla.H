#ifndef ANABeamPla_H
#define ANABeamPla_H

#include "AnaModule.H"

#include "TArtBigRIPSParameters.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibIC.hh"

class AnaBeamPla : public AnaModule {

public:

  AnaBeamPla();
  virtual ~AnaBeamPla();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
  
  TArtBigRIPSParameters* fBeamPlaParameters;
  
  TArtCalibPlastic* fCalibPlastic;
  TArtCalibIC* fCalibIC;

  Double_t dTRaw2X(Double_t dTRaw);
  Double_t DTtoX(Double_t dT);

  // for trees

  Double_t f3PlaTLRaw, f3PlaTRRaw;
  Double_t f5PlaTLRaw, f5PlaTRRaw;
  Double_t f7PlaTLRaw, f7PlaTRRaw;
  Double_t f13Pla1TLRaw, f13Pla1TRRaw;
  Double_t f13Pla2TLRaw, f13Pla2TRRaw;
  Double_t f3PlaTL, f3PlaTR; // Cal
  Double_t f3PlaQL, f3PlaQR; // Raw
  Double_t f5PlaTL, f5PlaTR; // Cal
  Double_t f5PlaQL, f5PlaQR; // Raw
  Double_t f7PlaTL, f7PlaTR; // Cal
  Double_t f7PlaQL, f7PlaQR; // Raw
  Double_t f13Pla1TL, f13Pla1TR; // Cal
  Double_t f13Pla1QL, f13Pla1QR; // Raw
  Double_t f13Pla2TL, f13Pla2TR; // Cal
  Double_t f13Pla2QL, f13Pla2QR; // Raw
  Double_t SBV1TL, SBV1TR; // Cal
  Double_t SBV1QL, SBV1QR; // Cal
  Double_t SBV2TL, SBV2TR; // Cal
  Double_t SBV2QL, SBV2QR; // Cal

  Double_t f5PlaTLSlew, f5PlaTRSlew;

  
  Double_t f5X; //

  Double_t icbEloss; //

};

#endif
