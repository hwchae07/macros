#include "AnaNeuLAND.H"

#include "TArtNEBULAPla.hh"
#include "TArtNeuLANDPla.hh"

#include "TTree.h"

AnaNeuLAND::AnaNeuLAND():
  AnaModule("NeuLAND"),
  fNeuLANDParameters(NULL),
  fCalibNeuLAND(NULL),
  fCalibNeuLANDVETO(NULL) {
  neulandNum = 0;
  neulandID = new Int_t[400];
  neulandQU = new Double_t[400];
  neulandQD = new Double_t[400];
  neulandTU = new Double_t[400];
  neulandTD = new Double_t[400];
  neulandTA = new Double_t[400];
  neulandX = new Double_t[400];
  neulandY = new Double_t[400];
  neulandZ = new Double_t[400];

  neulandTURaw = new Double_t[400];
  neulandTDRaw = new Double_t[400];
  //tofOffset = new Double_t[400];
}

AnaNeuLAND::~AnaNeuLAND(){    
  DeleteAll();
  delete neulandID;
  delete neulandQU;
  delete neulandQD;
  delete neulandTU;
  delete neulandTD;
  delete neulandTA;
  delete neulandX;
  delete neulandY;
  delete neulandZ;

  delete neulandTURaw;
  delete neulandTDRaw;
  //delete tofOffset;
}

void AnaNeuLAND::InitParameter(){
  fNeuLANDParameters = TArtSAMURAIParameters::Instance();
  fNeuLANDParameters->LoadParameter((char*)"db/NEULAND.xml");
  fNeuLANDParameters->LoadParameter((char*)"db/NEULANDVETO.xml");
  /*
  std::ifstream tof("data/neuland_tof_offset.dat");
  Double_t temp;
  for (Int_t i = 0 ; i < 400 ; i++)
    tof >> tofOffset[i];
  tof.close();
  */
  parLoaded = true;}

void AnaNeuLAND::InitDetector(){
  if (!parLoaded) return;
  fCalibNeuLAND = new TArtCalibNeuLAND;
  fCalibNeuLANDVETO = new TArtCalibNeuLANDVETO;
  detLoaded = true;}

void AnaNeuLAND::Analysis(){
  if (!detLoaded) return;
  anaFlag = true;

  fCalibNeuLAND->ClearData();
  fCalibNeuLANDVETO->ClearData();

  fCalibNeuLAND->ReconstructData();
  fCalibNeuLANDVETO->ReconstructData();

  std::ifstream fin("../../dat/neuland_gamma_offset.dat");
  Double_t neuland_offset;
  fin>>neuland_offset;
  fin.close();
  
  neulandNum = 0;  
  for (Int_t i = 0 ; i < 9 ; i++)
    neulandMult[i] = 0;
  for (Int_t i = 0 ; i < fCalibNeuLAND->GetNumNeuLANDPla() ; i++){
    TArtNeuLANDPla* pla = fCalibNeuLAND->GetNeuLANDPla(i);
    if (pla && pla->GetBothFired()){
      neulandID[neulandNum] = pla->GetID();
      neulandMult[int((neulandID[neulandNum]-1)/50)+1]++;
      neulandQU[neulandNum] = pla->GetQCal(0);
      neulandQD[neulandNum] = pla->GetQCal(1);
      neulandTU[neulandNum] = pla->GetTCal(0) - neuland_offset;// - tofOffset[neulandID[neulandNum]-1];
      neulandTD[neulandNum] = pla->GetTCal(1) - neuland_offset;// - tofOffset[neulandID[neulandNum]-1];
      //std::cout<<neuland_offset<<std::endl;
      neulandTA[neulandNum] = (neulandTU[neulandNum]+neulandTD[neulandNum])/2.;
      neulandX[neulandNum] = pla->GetX();
      neulandY[neulandNum] = pla->GetY();
      //std::cout << pla->GetZ() << std::endl;
      neulandZ[neulandNum] = pla->GetZ();
      neulandNum++;
    }
  }

  for (Int_t i = 0 ; i < fCalibNeuLANDVETO->GetNumNeuLANDVETOPla() ; i++){
    TArtNEBULAPla* pla = fCalibNeuLANDVETO->GetNeuLANDVETOPla(i);
    if (pla->GetHit()){
      if (pla->GetID() == 9) continue;
      neulandID[neulandNum] = pla->GetID()+400;
      neulandMult[0]++;
      neulandQU[neulandNum] = pla->GetQUCal();
      neulandQD[neulandNum] = pla->GetQDCal();
      neulandTU[neulandNum] = pla->GetTUSlw();
      neulandTD[neulandNum] = pla->GetTDSlw();
      neulandTA[neulandNum] = (pla->GetTUSlw()+pla->GetTDSlw())/2.;
      neulandX[neulandNum] = pla->GetPos(0);
      neulandY[neulandNum] = pla->GetPos(1);
      neulandZ[neulandNum] = pla->GetPos(2);
      neulandNum++;
    }
  }
  if (neulandNum < 1)
    anaFlag = false; 
}

void AnaNeuLAND::DeleteAll(){
  if (fCalibNeuLAND)    delete fCalibNeuLAND;
  if (fCalibNeuLANDVETO)    delete fCalibNeuLANDVETO;
}

void AnaNeuLAND::SetTree(){ 
  AnaModule::SetTree();
  
  tree->Branch("neulandNum",&neulandNum,"neulandNum/I");
  tree->Branch("neulandMult",neulandMult,"neulandMult[9]/I");
  tree->Branch("neulandID",neulandID,"neulandID[neulandNum]/I");

  tree->Branch("neulandQU",neulandQU,"neulandQU[neulandNum]/D");
  tree->Branch("neulandQD",neulandQD,"neulandQD[neulandNum]/D");
  tree->Branch("neulandTU",neulandTU,"neulandTU[neulandNum]/D");
  tree->Branch("neulandTD",neulandTD,"neulandTD[neulandNum]/D");
  tree->Branch("neulandTA",neulandTA,"neulandTA[neulandNum]/D");
  tree->Branch("neulandX",neulandX,"neulandX[neulandNum]/D");
  tree->Branch("neulandY",neulandY,"neulandY[neulandNum]/D");
  tree->Branch("neulandZ",neulandZ,"neulandZ[neulandNum]/D");

  tree->Branch("neulandTURaw",neulandTURaw,"neulandTURaw[neulandNum]/D");
  tree->Branch("neulandTDRaw",neulandTDRaw,"neulandTDRaw[neulandNum]/D");
}








