#include "AnaBeamPla.H"

#include "TTree.h"

#include "TArtPlastic.hh"
#include "TArtIC.hh"

AnaBeamPla::AnaBeamPla():
  AnaModule("BeamPla"),
  fBeamPlaParameters(NULL),
  fCalibPlastic(NULL),
  fCalibIC(NULL){
  ;}

AnaBeamPla::~AnaBeamPla(){    
  DeleteAll();}

void AnaBeamPla::InitParameter(){
  fBeamPlaParameters = TArtBigRIPSParameters::Instance();
  fBeamPlaParameters->LoadParameter((char*)"db/SAMURAIPlastic.xml");
  fBeamPlaParameters->LoadParameter((char*)"db/BigRIPSIC.xml");
  parLoaded = true;}

void AnaBeamPla::InitDetector(){
  if (!parLoaded) return;

  fCalibPlastic = new TArtCalibPlastic;
  fCalibIC = new TArtCalibIC;
  detLoaded = true;}

void AnaBeamPla::Analysis(){
  anaFlag = true;
  if (!detLoaded) return;

  fCalibPlastic->ClearData();
  fCalibPlastic->ReconstructData();

  fCalibIC->ClearData();
  fCalibIC->ReconstructData();

  TArtPlastic* pla;
  pla = fCalibPlastic->FindPlastic((char*)"F3pl");
  if (pla)
    {
      if(pla->GetTLRaw() >0 &&
	 pla->GetTRRaw() >0 &&
	 pla->GetQLRaw() > 0 &&
	 pla->GetQRRaw() > 0 &&
	 pla->GetTimeL() > -999 &&
	 pla->GetTimeR() > -999)
	{
	  f3PlaTLRaw = pla->GetTLRaw();
	  f3PlaTRRaw = pla->GetTRRaw();
	  f3PlaTL = pla->GetTimeL();
	  f3PlaTR = pla->GetTimeR();
	  f3PlaQL = pla->GetQLRaw();
	  f3PlaQR = pla->GetQRRaw();
	}
    else anaFlag = false; 
    // dT3 cut
    //    if (f3PlaTL - f3PlaTR < -5 ||
    //    	f3PlaTL - f3PlaTR > 5)
    //          anaFlag = false;
    } 
  else 
    anaFlag = false;

  pla = fCalibPlastic->FindPlastic((char*)"F5pl");
  if (pla) {
    if(pla->GetTLRaw() >0 &&
       pla->GetTRRaw() >0 &&
       pla->GetQLRaw() > 0 &&
       pla->GetQRRaw() > 0 &&
       pla->GetTimeL() > -999 &&
       pla->GetTimeR() > -999){
      f5PlaTLRaw = pla->GetTLRaw();
      f5PlaTRRaw = pla->GetTRRaw();
      f5PlaTLSlew = pla->GetTimeLSlew();
      f5PlaTRSlew = pla->GetTimeRSlew();      
      f5PlaTL = pla->GetTimeL();
      f5PlaTR = pla->GetTimeR();
      f5PlaQL = pla->GetQLRaw();
      f5PlaQR = pla->GetQRRaw();}
    else anaFlag = false;
    f5X = dTRaw2X(f5PlaTLRaw - f5PlaTRRaw);
    //f5X = DTtoX(f5PlaTL-f5PlaTR);
    // dT5 cut
    //    if (f5PlaTL - f5PlaTR < -20 ||
    //    	f5PlaTL - f5PlaTR > 20)
    //      anaFlag = false;
  }
  else
    anaFlag = false;
  
  pla = fCalibPlastic->FindPlastic((char*)"F7pl");
  if (pla) {
    if(pla->GetTLRaw() > 0 &&
       pla->GetTRRaw() > 0 &&
       pla->GetQLRaw() > 0 &&
       pla->GetQRRaw() > 0 &&
       pla->GetTimeL() > -999 &&
       pla->GetTimeR() > -999){
      f7PlaTLRaw = pla->GetTLRaw();
      f7PlaTRRaw = pla->GetTRRaw();
      f7PlaTL = pla->GetTimeL();
      f7PlaTR = pla->GetTimeR();
      f7PlaQL = pla->GetQLRaw();
      f7PlaQR = pla->GetQRRaw();}
    else anaFlag = false;
    // dT7 cut
    //    if (f7PlaTL - f7PlaTR < -15 ||
    //	f7PlaTL - f7PlaTR > 15)
    //      anaFlag = false;
  }
  else
    anaFlag = false;

  pla = fCalibPlastic->FindPlastic((char*)"F13pl-1");
  if (pla) {
    if(pla->GetTLRaw() > 0 &&
       pla->GetTRRaw() > 0 &&
       pla->GetQLRaw() > 0 &&
       pla->GetQRRaw() > 0 &&
       pla->GetTimeL() > -999 &&
       pla->GetTimeR() > -999){
      f13Pla1TLRaw = pla->GetTLRaw();
      f13Pla1TRRaw = pla->GetTRRaw();
      f13Pla1TL = pla->GetTimeL();
      f13Pla1TR = pla->GetTimeR();
      f13Pla1QL = pla->GetQLRaw();
      f13Pla1QR = pla->GetQRRaw();}
    else anaFlag = false; }
  else
    anaFlag = false;

  pla = fCalibPlastic->FindPlastic((char*)"F13pl-2");
  if (pla) {
    if(pla->GetTLRaw() > 0 &&
       pla->GetTRRaw() > 0 &&
       pla->GetQLRaw() > 0 &&
       pla->GetQRRaw() > 0 &&
       pla->GetTimeL() > -999 &&
       pla->GetTimeR() > -999){
      f13Pla2TLRaw = pla->GetTLRaw();
      f13Pla2TRRaw = pla->GetTRRaw();
      f13Pla2TL = pla->GetTimeL();
      f13Pla2TR = pla->GetTimeR();
      f13Pla2QL = pla->GetQLRaw();
      f13Pla2QR = pla->GetQRRaw();}
    else anaFlag = false; }
  else
    anaFlag = false;

  //ICB
  TArtIC *ic = fCalibIC->FindIC((char*)"ICB");
  
  if (ic)
    icbEloss = ic->GetEnergyAvSum();	
  else
    anaFlag = false;

  
}

void AnaBeamPla::DeleteAll(){
  if (fCalibPlastic)    delete fCalibPlastic;
  if (fCalibIC)    delete fCalibIC;
}

void AnaBeamPla::SetTree(){ 
  AnaModule::SetTree();
  tree->Branch("f3PlaTLRaw",&f3PlaTLRaw,"f3PlaTLRaw/D");
  tree->Branch("f3PlaTRRaw",&f3PlaTRRaw,"f3PlaTRRaw/D");
  tree->Branch("f5PlaTLRaw",&f5PlaTLRaw,"f5PlaTLRaw/D");
  tree->Branch("f5PlaTRRaw",&f5PlaTRRaw,"f5PlaTRRaw/D");
  tree->Branch("f7PlaTLRaw",&f7PlaTLRaw,"f7PlaTLRaw/D");
  tree->Branch("f7PlaTRRaw",&f7PlaTRRaw,"f7PlaTRRaw/D");
  tree->Branch("f13Pla1TLRaw",&f13Pla1TLRaw,"f13Pla1TLRaw/D");
  tree->Branch("f13Pla1TRRaw",&f13Pla1TRRaw,"f13Pla1TRRaw/D");
  tree->Branch("f13Pla2TLRaw",&f13Pla2TLRaw,"f13Pla2TLRaw/D");
  tree->Branch("f13Pla2TRRaw",&f13Pla2TRRaw,"f13Pla2TRRaw/D");
  
  tree->Branch("f3PlaTL",&f3PlaTL,"f3PlaTL/D");
  tree->Branch("f3PlaTR",&f3PlaTR,"f3PlaTR/D");
  tree->Branch("f3PlaQL",&f3PlaQL,"f3PlaQL/D");
  tree->Branch("f3PlaQR",&f3PlaQR,"f3PlaQR/D");
  tree->Branch("f5PlaTL",&f5PlaTL,"f5PlaTL/D");
  tree->Branch("f5PlaTR",&f5PlaTR,"f5PlaTR/D");
  tree->Branch("f5PlaQL",&f5PlaQL,"f5PlaQL/D");
  tree->Branch("f5PlaQR",&f5PlaQR,"f5PlaQR/D");
  tree->Branch("f7PlaTL",&f7PlaTL,"f7PlaTL/D");
  tree->Branch("f7PlaTR",&f7PlaTR,"f7PlaTR/D");
  tree->Branch("f7PlaQL",&f7PlaQL,"f7PlaQL/D");
  tree->Branch("f7PlaQR",&f7PlaQR,"f7PlaQR/D");
  tree->Branch("f13Pla1TL",&f13Pla1TL,"f13Pla1TL/D");
  tree->Branch("f13Pla1TR",&f13Pla1TR,"f13Pla1TR/D");
  tree->Branch("f13Pla1QL",&f13Pla1QL,"f13Pla1QL/D");
  tree->Branch("f13Pla1QR",&f13Pla1QR,"f13Pla1QR/D");
  tree->Branch("f13Pla2TL",&f13Pla2TL,"f13Pla2TL/D");
  tree->Branch("f13Pla2TR",&f13Pla2TR,"f13Pla2TR/D");
  tree->Branch("f13Pla2QL",&f13Pla2QL,"f13Pla2QL/D");
  tree->Branch("f13Pla2QR",&f13Pla2QR,"f13Pla2QR/D");
  
  tree->Branch("f5PlaTLSlew",&f5PlaTLSlew,"f5PlaTLSlew/D");
  tree->Branch("f5PlaTRSlew",&f5PlaTRSlew,"f5PlaTRSlew/D");
  
  tree->Branch("f5X",&f5X,"f5X/D");
  tree->Branch("icbEloss",&icbEloss,"icbEloss/D");

  
  tree->SetAlias("T3","(f3PlaTL+f3PlaTR)/2");
  tree->SetAlias("T7","(f7PlaTL+f7PlaTR)/2");
  tree->SetAlias("T13","(f13Pla1TL+f13Pla1TR+f13Pla2TL+f13Pla2TR)/4");
  tree->SetAlias("TOF3_13","T13-T3");
  tree->SetAlias("TOF7_13","T13-T7");
  tree->SetAlias("TOF3_7","T7-T3");
  tree->SetAlias("dT3","f3PlaTL-f3PlaTR");
  tree->SetAlias("dT5","f5PlaTL-f5PlaTR");
  tree->SetAlias("dT5Raw","f5PlaTLRaw-f5PlaTRRaw");
  tree->SetAlias("dT5Slew","f5PlaTLSlew-f5PlaTRSlew");
  tree->SetAlias("dT7","f7PlaTL-f7PlaTR");
  tree->SetAlias("dT13","f13Pla1TL-f13Pla1TR");
  tree->SetAlias("Q13","sqrt(f13Pla1QL*f13Pla1QR)");
  tree->SetAlias("Q3","sqrt(f3PlaQL*f3PlaQR)");
  tree->SetAlias("Q7","sqrt(f7PlaQL*f7PlaQR)");
  tree->SetAlias("Q5","sqrt(f5PlaQL*f5PlaQR)");
}

Double_t AnaBeamPla::dTRaw2X(Double_t dTRaw)
{

  Double_t par1[6] = {-104.582,
		      -0.189221,
		      0.00166446,
		      -1.93474e-06,
		      -1.66557e-08,
		      -1.62186e-11};
  Double_t par2[6] = {-102.967,
		      -0.179801,
		      0.0019141,
		      -8.03565e-08,
		      -1.17903e-08,
		      -1.21168e-11};
  Double_t par3[6] = {-107.217,
		      -0.165188,
		      0.00259371,
		      5.11836e-06,
		      2.54618e-09,
		      1.25692e-12};
  Double_t par4[6] = {-107.417,
		      -0.175539,
		      0.00253763,
		      5.06046e-06,
		      2.83092e-09,
		      1.8012e-12};
  Double_t x1 = 0;
  for(Int_t i=0 ; i<6 ; i++)
    x1+=par1[i] * TMath::Power(dTRaw,i);
  Double_t x2 = 0;
  for(Int_t i=0 ; i<6 ; i++)
    x2+=par2[i] * TMath::Power(dTRaw,i);
  Double_t x3 = 0;
  for(Int_t i=0 ; i<6 ; i++)
    x3+=par3[i] * TMath::Power(dTRaw,i);
  Double_t x4 = 0;
  for(Int_t i=0 ; i<6 ; i++)
    x4+=par4[i] * TMath::Power(dTRaw,i);

  Double_t X1 = (x1 + x2) / 2.;
  Double_t X2 = (x3 + x4) / 2.;

  Double_t X = (450.*X2 + 200.*X1)/(450.+250.);

  return X;
}


Double_t AnaBeamPla::DTtoX(Double_t dT){

  Double_t par[7] = {
    112.553,
    -39.0031,
    -14.2877,
    5.74416,
   -1.08954,
    0.0811626,
    0};
  Double_t x = 0;
  Double_t temp = 1;
  for (Int_t i = 0 ; i < 7 ; i++){
    x += (par[i]*temp);
    temp *= dT; }
  return x;
}
