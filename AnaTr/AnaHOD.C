#include <fstream>
#include <iostream>

#include "AnaHOD.H"

#include "TTree.h"

#include "TArtHODPla.hh"

#include "TMath.h"
#include "TString.h"

AnaHOD::AnaHOD():
  AnaModule("HOD"),
  fHODParameters(NULL),
  fCalibHODPla(NULL),
  hodNum(0) {
  hodID = new Int_t[24];
  hodPlaQU = new Double_t[24];
  hodPlaQD = new Double_t[24];
  hodPlaTU = new Double_t[24];
  hodPlaTD = new Double_t[24];
  hodPlaTA = new Double_t[24];
  hodTURaw = new Double_t[24];
  hodTDRaw = new Double_t[24];

  hodQURaw = new Double_t[24];
  hodQDRaw = new Double_t[24];
  
  hodQ = new Double_t[24];
}

AnaHOD::~AnaHOD(){    
  DeleteAll();
  delete hodID;
  delete hodPlaQU;
  delete hodPlaQD;
  delete hodPlaTU;
  delete hodPlaTD;
  delete hodPlaTA;
  delete hodTURaw;
  delete hodTDRaw;

  delete hodQURaw;
  delete hodQDRaw;
  
  delete hodQ;
}

void AnaHOD::InitParameter(){
  fHODParameters = TArtSAMURAIParameters::Instance();
  fHODParameters->LoadParameter((char*)"db/SAMURAIHOD.xml");

  parLoaded = true;}

void AnaHOD::InitDetector(){
  if (!parLoaded) return;
  fCalibHODPla = new TArtCalibHODPla;
  detLoaded = true;}

void AnaHOD::Analysis(){
  anaFlag = true;
  if (!detLoaded) return;

  fCalibHODPla->ClearData();
  fCalibHODPla->ReconstructData();

  std::ifstream fin1("../../dat/hod_offset_brho_scan.dat");
  Int_t dummy[25];
  Double_t hodT_offset[25];
  for(Int_t i=0;i<=24;i++)
    fin1>>dummy[i]>>hodT_offset[i];
  fin1.close();
  
  std::ifstream fin2("../../dat/physics_run_hodq.dat");
  Double_t hodQ_offset[25];
  for(Int_t id = 1;id<=24;id++)
    fin2>>hodQ_offset[id];
  fin2.close();  

  //  hodNum = fCalibHODPla->GetNumHODPla();
  //  if ( hodNum < 1 ) {
  //    anaFlag = false; return; }
  hodNum = 0;
  for (Int_t i = 0 ; i < fCalibHODPla->GetNumHODPla() ; i++){
    TArtHODPla* pla = fCalibHODPla->GetHODPla(i);
    hodID[hodNum] = pla->GetID();
    hodPlaQU[hodNum] = pla->GetQUCal();
    hodPlaQD[hodNum] = pla->GetQDCal();
    hodPlaTU[hodNum] = pla->GetTimeUSlew();// + parTOF[0] + parTOF[pla->GetID()];
    hodPlaTD[hodNum] = pla->GetTimeDSlew();// + parTOF[0] + parTOF[pla->GetID()];
    hodPlaTA[hodNum] = pla->GetTimeSlew() + hodT_offset[0] + hodT_offset[pla->GetID()];
    hodTURaw[hodNum] = pla->GetTURaw();
    hodTDRaw[hodNum] = pla->GetTDRaw();

    hodQURaw[hodNum] = pla->GetQURaw();
    hodQDRaw[hodNum] = pla->GetQDRaw();
    

    hodQ[hodNum] = TMath::Sqrt(hodPlaQU[hodNum] * hodPlaQD[hodNum]) * hodQ_offset[hodID[hodNum]];
    
    hodNum++;}
  


  /*  
  //2017.02.16 Hyunwoo//
  //ID sort//
  Int_t index[24]={0,};
  TMath::Sort(hodNum,hodID,index,false);
  for(Int_t i=0 ; i<hodNum ; i++)
    hodID[i] = hodID[index[i]];
  //2017.02.16 Hyunwoo//
  */

  
  if (!hodNum) anaFlag = false;
}

void AnaHOD::DeleteAll(){
  if (fCalibHODPla)    delete fCalibHODPla;
}

void AnaHOD::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("hodNum",&hodNum,"hodNum/I");
  tree->Branch("hodID",hodID,"hodID[hodNum]/I");  
  tree->Branch("hodQU",hodPlaQU,"hodQU[hodNum]/D");
  tree->Branch("hodQD",hodPlaQD,"hodQD[hodNum]/D");
  tree->Branch("hodTU",hodPlaTU,"hodTU[hodNum]/D");
  tree->Branch("hodTD",hodPlaTD,"hodTD[hodNum]/D");
  tree->Branch("hodTA",hodPlaTA,"hodTA[hodNum]/D");
  tree->Branch("hodTURaw",hodTURaw,"hodTURaw[hodNum]/D");
  tree->Branch("hodTDRaw",hodTDRaw,"hodTDRaw[hodNum]/D");
  
  tree->Branch("hodQURaw",hodQURaw,"hodQURaw[hodNum]/D");  
  tree->Branch("hodQDRaw",hodQDRaw,"hodQDRaw[hodNum]/D");


  tree->Branch("hodQ",hodQ,"hodQ[hodNum]/D");
  
  tree->SetAlias("hodQA","sqrt(hodQU*hodQD)");
  

}

