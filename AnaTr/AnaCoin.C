#include "AnaCoin.H"

#include "TTree.h"
#include "TClonesArray.h"

#include "TArtPlastic.hh"
#include "TArtStoreManager.hh"
#include "TArtEventInfo.hh"

AnaCoin::AnaCoin():
  AnaModule("Coin"),
  fCalibCoin(NULL) {
  ;}

AnaCoin::~AnaCoin(){    
  DeleteAll();}

void AnaCoin::InitParameter(){
  parLoaded = true;}

void AnaCoin::InitDetector(){
  if (!parLoaded) return;
  fCalibCoin = new TArtCalibCoin;
  detLoaded = true;}

void AnaCoin::Analysis(){
  anaFlag = true;
  if (!detLoaded) return;

  fCalibCoin->ClearData();
  fCalibCoin->LoadData();

  TArtStoreManager* fStoreManager;
  fStoreManager = TArtStoreManager::Instance();
  TArtEventInfo* eventInfo = (TArtEventInfo*)((TClonesArray*)fStoreManager->FindDataContainer("EventInfo"))->At(0);
  if (eventInfo) coinTrigger = eventInfo->GetTriggerBit() & 0xFF;
  else anaFlag = false;
}

void AnaCoin::DeleteAll(){
  if (fCalibCoin)    delete fCalibCoin;
}

void AnaCoin::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("coinTrigger",&coinTrigger,"coinTrigger/I");
}

