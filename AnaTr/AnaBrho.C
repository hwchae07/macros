#include <iostream>

//#include "SAMURAIMatrix.cc"

#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TMath.h"

#include "Julien.H"

using namespace std;

void AnaBrho(Int_t runNum) {

  Julien* gib = new Julien();
  Bool_t goodFlag = gib->ReadParameterMultiDimFit();

  if (!goodFlag) {
    cout << "No way!" << endl;
    return;}

  TFile* file = new TFile(Form("../../root/run%04d.root.FDC",runNum),"READ");
  TTree* tree;
  file->GetObject("FDC",tree);
  tree->BuildIndex("EventNum");
  tree->AddFriend("BDC",Form("../../root/run%04d.root.BDC",runNum));
  tree->GetFriend("BDC")->BuildIndex("EventNum");

  Int_t RunNum, EventNum;
  Double_t fdc1X, fdc1Y, fdc1A, fdc1B;
  Double_t fdc2X, fdc2Y, fdc2A, fdc2B;
  Double_t bdc1X, bdc1Y, bdc2X, bdc2Y;

  tree->SetBranchAddress("RunNum",&RunNum);
  tree->SetBranchAddress("EventNum",&EventNum);
  tree->SetBranchAddress("fdc1X",&fdc1X);
  tree->SetBranchAddress("fdc1Y",&fdc1Y);
  tree->SetBranchAddress("fdc1A",&fdc1A);
  tree->SetBranchAddress("fdc1B",&fdc1B);
  tree->SetBranchAddress("fdc2X",&fdc2X);
  tree->SetBranchAddress("fdc2Y",&fdc2Y);
  tree->SetBranchAddress("fdc2A",&fdc2A);
  tree->SetBranchAddress("fdc2B",&fdc2B);
  tree->SetBranchAddress("bdc1X",&bdc1X);
  tree->SetBranchAddress("bdc1Y",&bdc1Y);
  tree->SetBranchAddress("bdc2X",&bdc2X);
  tree->SetBranchAddress("bdc2Y",&bdc2Y);

  Double_t brho, fl;
  TFile *file2 = new TFile(Form("../../root/run%04d.root.Brho",runNum),"RECREATE");
  TTree *tree2 = new TTree("Brho","Brho");
  tree2->Branch("RunNum",&RunNum,"RunNum/I");
  tree2->Branch("EventNum",&EventNum,"EventNum/I");
  tree2->Branch("brho",&brho,"brho/D");
  tree2->Branch("fl",&fl,"fl/D");

  for (Int_t i = 0 ; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if (tree->GetEntryWithIndex(EventNum) < 0) continue;
    
    Double_t tgtX = bdc1X + 1853. / 1000. * (bdc2X - bdc1X);
    Double_t tgtY = bdc1Y + 1853. / 1000. * (bdc2Y - bdc1Y);

    Double_t in[6];
    in[0] = fdc1X;
    in[1] = -fdc1Y;
    in[2] = TMath::ATan2( fdc1X - tgtX, 1680.);
    in[3] = TMath::ATan2(-fdc1Y - tgtY, 1680.);
    in[4] = fdc2X;
    in[5] = TMath::ATan(fdc2A);
    
    fl   = gib->GetPathLengthFitResult(in);
    brho = gib->GetRigidityFitResult(in);
    tree2->Fill();
  }
  tree2->Write();
  file2->Close();
}
