#include "AnaNEBslew.H"

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCut.h"

void AnaNEBslew(Int_t runNum)
{
  TFile *file = new TFile(Form("../../root/run%04d.root.BeamPID",runNum),"READ");
  TTree *tree;
  file->GetObject("BeamPID",tree);
  tree->BuildIndex("EventNum");
  tree->AddFriend("NEBULA",Form("../../root/run%04d.root.NEBULA",runNum));
  tree->GetFriend("NEBULA")->BuildIndex("EventNum");

  Int_t RunNum,EventNum;
  Double_t tofc713;
  Double_t beamAoZ, beamZ;

  Int_t nebulaNum;
  Int_t* nebulaID;
  Double_t* nebulaTU;
  Double_t* nebulaTD;
  Double_t* nebulaTA;
  Double_t* nebulaQUPed;
  Double_t* nebulaQDPed;

  nebulaID = new Int_t[144];
  nebulaTU = new Double_t[144];
  nebulaTD = new Double_t[144];
  nebulaTA = new Double_t[144];
  nebulaQUPed = new Double_t[144];
  nebulaQDPed = new Double_t[144];

  tree->SetBranchAddress("RunNum",&RunNum);
  tree->SetBranchAddress("EventNum",&EventNum);

  tree->SetBranchAddress("tofc713",&tofc713);  

  tree->SetBranchAddress("beamAoZ",&beamAoZ);
  tree->SetBranchAddress("beamZ",&beamZ);
  
  tree->SetBranchAddress("nebulaNum",&nebulaNum);
  tree->SetBranchAddress("nebulaID",nebulaID);
  tree->SetBranchAddress("nebulaTU",nebulaTU);
  tree->SetBranchAddress("nebulaTD",nebulaTD);
  tree->SetBranchAddress("nebulaTA",nebulaTA);
  tree->SetBranchAddress("nebulaQUPed",nebulaQUPed);
  tree->SetBranchAddress("nebulaQDPed",nebulaQDPed);

  
  TFile *file2 = new TFile(Form("../../root/run%04d.root.NEBslew",runNum),"RECREATE");
  TTree *tree2 = new TTree("NEBslew","NEBslew");
  tree2->Branch("RunNum",&RunNum,"RunNum/I");
  tree2->Branch("EventNum",&EventNum,"EventNum/I");
  tree2->Branch("tofc713",&tofc713,"tofc713/D");
  tree2->Branch("beamAoZ",&beamAoZ,"beamAoZ/D");
  tree2->Branch("beamZ",&beamZ,"beamZ/D");

  tree2->Branch("nebulaNum",&nebulaNum,"nebulaNum/I");
  tree2->Branch("nebulaID",nebulaID,"nebulaID[nebulaNum]/I");
  tree2->Branch("nebulaTU",nebulaTU,"nebulaTU[nebulaNum]/D");
  tree2->Branch("nebulaTD",nebulaTD,"nebulaTD[nebulaNum]/D");
  tree2->Branch("nebulaTA",nebulaTA,"nebulaTA[nebulaNum]/D");
  tree2->Branch("nebulaQUPed",nebulaQUPed,"nebulaQUPed[nebulaNum]/D");
  tree2->Branch("nebulaQDPed",nebulaQDPed,"nebulaQDPed[nebulaNum]/D");
  
  for(Int_t i=0;i<tree->GetEntries();i++)
    {
      if(!(i%1000))
        {
          std::cout<<"Progress rate(%): " << (Double_t)i/(Double_t)tree->GetEntries()*100 << "\r";
          std::cout.flush();
        }

      tree->GetEntry(i);
      if(tree->GetEntryWithIndex(EventNum)<0) continue;

      /*
	for(Int_t j=0;j<nebulaNum;j++)
	{
	nebulaID[j] = nebulaID[j];
	nebulaTU[j] = nebulaTU[j];
	nebulaTD[j] = nebulaTD[j];
	nebulaTA[j] = nebulaTA[j];
	nebulaQUPed[j] = nebulaQUPed[j];
	nebulaQDPed[j] = nebulaQDPed[j];
	}
      */
      tree2->Fill();
    }

  delete nebulaID;
  delete nebulaTU;
  delete nebulaTD;
  delete nebulaTA;
  delete nebulaQUPed;
  delete nebulaQDPed;

  /*
    delete neutID;
    delete neutTA;
    delete neutX;
    delete neutY;
    delete neutZ;
    delete neutIndex;
    delete neutNEB;
    delete neutVETO;
    delete neutFL1;
    delete neutTOF1;
    delete neutBeta1;
  */

  tree2->Print();
  tree2->Write();
  file2->Close();
}