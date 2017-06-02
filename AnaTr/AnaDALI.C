#include <iostream>
#include <fstream>

#include "AnaDALI.H"

#include "TArtDALINaI.hh"

#include "TTree.h"

AnaDALI::AnaDALI():
  AnaModule("DALI"),
  fDALIParameters(NULL),
  fCalibDALI(NULL) {
  daliNum = 0;
  daliID = new Int_t[140];
  daliTheta = new Double_t[140];
  daliE = new Double_t[140];
  daliT = new Double_t[140];
}

AnaDALI::~AnaDALI(){    
  DeleteAll();
  delete daliID;
  delete daliTheta;
  delete daliE;
  delete daliT;
}

void AnaDALI::InitParameter(){
  fDALIParameters = TArtDALIParameters::Instance();
  fDALIParameters->LoadParameter((char*)"db/DALI.xml");

  parLoaded = true;}

void AnaDALI::InitDetector(){
  if (!parLoaded) return;
  fCalibDALI = new TArtCalibDALI;
  detLoaded = true;}

void AnaDALI::Analysis(){
  if (!detLoaded) return;
  anaFlag = true;

  fCalibDALI->ClearData();
  fCalibDALI->ReconstructData();

  daliNum = 0;

  for (Int_t i = 0 ; i < fCalibDALI->GetNumNaI() ; i++){
    TArtDALINaI* nai = fCalibDALI->GetNaI(i);
    if (nai && nai->GetID() < 141 && 
	nai->GetRawADC() > 0 && nai->GetRawADC() < 4095 &&
	nai->GetTimeOffseted() > 0){
      daliID[daliNum] = nai->GetID();
      daliTheta[daliNum] = nai->GetTheta();
      daliE[daliNum] = nai->GetEnergy();
      daliT[daliNum] = nai->GetTimeOffseted();
      daliNum++;
    }
  }
  if (daliNum < 1)
    anaFlag = false; 
}

void AnaDALI::DeleteAll(){
  if (fCalibDALI)    delete fCalibDALI;
}

void AnaDALI::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("daliNum",&daliNum,"daliNum/I");
  tree->Branch("daliID",daliID,"daliID[daliNum]/I");
  tree->Branch("daliTheta",daliTheta,"daliTheta[daliNum]/D");
  tree->Branch("daliE",daliE,"daliE[daliNum]/D");
  tree->Branch("daliT",daliT,"daliT[daliNum]/D");
}
