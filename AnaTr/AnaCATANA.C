#include <iostream>
#include <fstream>

#include "AnaCATANA.H"

#include "TArtCATANACsI.hh"

#include "TMath.h"
#include "TTree.h"

AnaCATANA::AnaCATANA():
  AnaModule("CATANA"),
  fCATANAParameters(NULL),
  fCalibCATANA(NULL) {
  catanaNum = 0;
  catanaID = new Int_t[100];
  catanaTheta = new Double_t[100];
  catanaE = new Double_t[100];
  catanaQ = new Double_t[100];
  catanaT = new Double_t[100];
  catanaRawTRef = new Double_t[100];
  catanaRawTDC = new Double_t[100];
}

AnaCATANA::~AnaCATANA(){    
  DeleteAll();
  delete catanaID;
  delete catanaTheta;
  delete catanaE;
  delete catanaQ;
  delete catanaT;
  delete catanaRawTDC;
  delete catanaRawTRef;
}

void AnaCATANA::InitParameter(){
  fCATANAParameters = TArtCATANAParameters::Instance();
  fCATANAParameters->LoadParameter((char*)"db/CATANA.xml");

  parLoaded = true;}

void AnaCATANA::InitDetector(){
  if (!parLoaded) return;
  fCalibCATANA = new TArtCalibCATANA;
  detLoaded = true;}

void AnaCATANA::Analysis(){
  if (!detLoaded) return;
  anaFlag = true;

  fCalibCATANA->ClearData();
  fCalibCATANA->ReconstructData();

  catanaNum = 0;

  for (Int_t i = 0 ; i < fCalibCATANA->GetCsIEntries() ; i++){
    TArtCATANACsI* csi = fCalibCATANA->GetCsI(i);
    if (csi && csi->GetID() < 141 && 
	csi->GetRawADC() > 0 && csi->GetRawADC() < 4095 &&
	csi->GetTimeOffseted() > 0){
      catanaID[catanaNum] = csi->GetID();
      //      catanaTheta[catanaNum] = csi->GetTheta();
      Double_t x =  csi->GetPositionX();
      Double_t y =  csi->GetPositionY();
      Double_t z =  csi->GetPositionZ();
      Double_t r = TMath::Sqrt(x*x+y*y);
      catanaTheta[catanaNum] = TMath::ATan(r/z);
      catanaE[catanaNum] = csi->GetEnergy();
      catanaQ[catanaNum] = csi->GetRawADC();
      catanaT[catanaNum] = csi->GetTimeOffseted();
      catanaRawTDC[catanaNum] = csi->GetRawTDC();
      catanaRawTRef[catanaNum] = csi->GetRawTRef();
      catanaNum++;
    }
  }
  if (catanaNum < 1)
    anaFlag = false; 
}

void AnaCATANA::DeleteAll(){
  if (fCalibCATANA)    delete fCalibCATANA;
}

void AnaCATANA::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("catanaNum",&catanaNum,"catanaNum/I");
  tree->Branch("catanaID",catanaID,"catanaID[catanaNum]/I");
  tree->Branch("catanaTheta",catanaTheta,"catanaTheta[catanaNum]/D");
  tree->Branch("catanaE",catanaE,"catanaE[catanaNum]/D");
  tree->Branch("catanaQ",catanaQ,"catanaQ[catanaNum]/D");
  tree->Branch("catanaT",catanaT,"catanaT[catanaNum]/D");
  tree->Branch("catanaRawTDC",catanaRawTDC,"catanaRawTDC[catanaNum]/D");
  tree->Branch("catanaRawTRef",catanaRawTRef,"catanaRawTRef[catanaNum]/D");
}
