#include "AnaTree.H"

#include <iostream>
#include <cstdlib>

#include "AnaModule.H"
//#include "AnaMINOS.H"
#include "AnaNeut.H"
#include "AnaNEBULA.H"
#include "AnaNeuLAND.H"
#include "AnaBeamPla.H"
#include "AnaBDC.H"
#include "AnaFDC.H"
#include "AnaHOD.H"
#include "AnaCoin.H"
#include "AnaDALI.H"
#include "AnaCATANA.H"
#include "AnaPPAC.H"
#include "AnaBrho.H"
#include "AnaSAMURAI.H"
#include "AnaBeamPID.H"
#include "AnaFragPID.H"
#include "AnaRelE.H"
#include "AnaNEBslew.H"
#include "AnaHPC.H"

#include "TArtStoreManager.hh"
#include "TArtEventStore.hh"

#include "TString.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TPRegexp.h"

using namespace std;

int main(int argc, char** argv){
  flagEnd = false;

  PrintTitle();
  // parsing
  if (argc < 2) {
    PrintUsage();
    return 0; }

  Int_t nMaxEvent = 0;
  TString ridfFileName(argv[1]);
  
  if (argc >= 3) nMaxEvent = atoi(argv[2]);

  // Load Modules
  AnaModule* modules[11] = { new AnaCoin(), new AnaBeamPla(), new AnaBDC(), new AnaHOD(), new AnaFDC(), new AnaNeut(), new AnaNEBULA(), new AnaNeuLAND(), new AnaCATANA(), new AnaPPAC(), new AnaHPC() };
  Int_t whichModule = SelectModule();
  if (!whichModule) return 0;

  /*
    if (whichModule == 11) { // for AnaBrho
    AnaBrho(atoi(argv[1]));
    return 0;}
  */
  if (whichModule == 12)
    {
      AnaSAMURAI(atoi(argv[1]));
      return 0;
    }

  if (whichModule == 13)
    {
      AnaBeamPID(atoi(argv[1]));
      return 0;
    }

  if (whichModule == 14)
    {
      AnaFragPID(atoi(argv[1]));
      return 0;
    }

  if (whichModule == 15)
    {
      AnaRelE(atoi(argv[1]));
      return 0;
    }

  if (whichModule == 16)
    {
      AnaNEBslew(atoi(argv[1]));
      return 0;
    }

  
  // Construction Part
  TArtStoreManager* fStoreManager;
  fStoreManager = TArtStoreManager::Instance();

  AnaModule* fAnaModule;
  fAnaModule = modules[whichModule-1];
  fAnaModule->InitParameter();
  fAnaModule->InitDetector();

  //  TObjArray* token = ridfFileName->Tokenize();
  TStringToken* token = new TStringToken(ridfFileName,"/");
  while (token->NextToken());
  token->Resize(token->Length()-5);
  //  const char* rootName = "test";
  TFile *rootFile = LoadRootFile((const char*)(*token),(const char*)fAnaModule->GetDetName(),false);
  fAnaModule->InitTree();
  fAnaModule->SetTree();
  
  // Calculate (Event by event)
  TArtEventStore* eventStore = new TArtEventStore;
  if (!eventStore->Open((const char*)ridfFileName)){
    // no ridf file
    exit(0); }

  Int_t nEvent = 0;

  cout << " ============================================= " << endl;
  cout << "              < Start analysis > " << endl;
  cout << " ============================================= " << endl;
  while (!flagEnd && eventStore->GetNextEvent() && (!nMaxEvent || nEvent < nMaxEvent)){
    signal(SIGINT, &SigHandler);

    if (!(nEvent % 1000)) {
      cout << "   Event number : " << nEvent << "\r";
      cout.flush();}
    fAnaModule->PreAnalysis();
    fAnaModule->Analysis();

    fAnaModule->FillTree();
    nEvent++;}

  cout << "   Event Number : " << nEvent << endl;
  cout << " ============================================= " << endl;
  cout << "               < End analysis > " << endl;
  cout << " ============================================= " << endl;
  cout << "              < Tree information > " << endl;
  fAnaModule->GetTree()->Print();
  cout << " ============================================= " << endl;
  cout << "               <Tree print done>               " <<endl;
  fAnaModule->GetTree()->Write();
  cout << "               <Tree write done>               " <<endl;
  cout << " ============================================= " << endl;
  rootFile->Close();
  delete rootFile;

  //delete fAnaModule;
  delete fStoreManager;
  return 0;
}


void PrintUsage(){

  cout << " USAGE :" << endl;
  cout << "  AnaTree [ridffile] ([# of events])" << endl;
  cout << "  by J. Hwang" << endl;}

void PrintTitle(){
  cout << " ============================================= " << endl;
  cout << "            AnaTree (11th Nov. 2016)" << endl;
  cout << " ============================================= " << endl;
}
int SelectModule(){
  cout << " ============================================= " << endl;
  cout << "              < Select a module > " << endl;
  cout << " ============================================= " << endl;
  cout << "   0. Exit" << endl;
  cout << "   1. Coin. Reg." << endl;
  cout << "   2. Beam Pla." << endl;
  cout << "   3. BDC" << endl;
  //  cout << "   4. MINOS" << endl;
  cout << "   4. HOD" << endl;
  cout << "   5. FDC" << endl;
  cout << "   6. Neut" << endl;
  cout << "   7. NEBULA" << endl;
  cout << "   8. NeuLAND" << endl;
  cout << "   9. CATANA" << endl;
  cout << "   10. PPAC" << endl;
  cout << "   11. NEBULA with HPC" <<endl;
  //  cout << "   11. Brho" << endl;
  cout << "   12. Brho (by Chae)" <<endl;
  cout << "   13. BeamPID" << endl;
  cout << "   14. FragPID" << endl;
  cout << "   15. RelE" <<endl;
  cout << "   16. NEBULA for slew" <<endl;
  cout << " ============================================= " << endl;
  Int_t n;

  do {
    cout << "  > ";  
    cin >> n;
    if (n >= 0 && n <= 16) return n;
  } while(1);

  return n;}

TFile* LoadRootFile(const char* name, const char* detName, bool readOnly){

  //TString rootName("root/");
  TString rootName("/home/hwchae07/Work/SAMURAI2016/root/");
  rootName += name;
  rootName += ".root.";
  rootName += detName;

  TFile* file;
  if (readOnly)
    file = new TFile((const char*)rootName,"READ");
  else
    file = new TFile((const char*)rootName,"RECREATE");

  return file;}
