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
  Double_t f13Pla1TL,f13Pla1TR,f13Pla2TL,f13Pla2TR;
  Double_t tofc713;
  Double_t beamAoZ, beamZ;

  Int_t nebulaNum;
  Int_t* nebulaID;
  Double_t* nebulaTURaw;
  Double_t* nebulaTDRaw;
  Double_t* nebulaTU;
  Double_t* nebulaTU1;
  Double_t* nebulaTU2;
  Double_t* nebulaTD;
  Double_t* nebulaTD1;
  Double_t* nebulaTD2;
  Double_t* nebulaTA;
  Double_t* nebulaTA1;
  Double_t* nebulaTA2;
  Double_t* nebulaQUPed;
  Double_t* nebulaQDPed;
  Double_t* nebulaQPed;
  Double_t* nebulaTOF;

  nebulaID = new Int_t[144];
  nebulaTURaw = new Double_t[144];
  nebulaTDRaw = new Double_t[144];
  nebulaTU = new Double_t[144];
  nebulaTU1 = new Double_t[144];
  nebulaTU2 = new Double_t[144];
  nebulaTD = new Double_t[144];
  nebulaTD1 = new Double_t[144];
  nebulaTD2 = new Double_t[144];
  nebulaTA = new Double_t[144];
  nebulaTA1 = new Double_t[144];
  nebulaTA2 = new Double_t[144];
  nebulaQUPed = new Double_t[144];
  nebulaQDPed = new Double_t[144];
  nebulaQPed = new Double_t[144];
  nebulaTOF = new Double_t[144];
  
  tree->SetBranchAddress("RunNum",&RunNum);
  tree->SetBranchAddress("EventNum",&EventNum);

  tree->SetBranchAddress("tofc713",&tofc713);  
  tree->SetBranchAddress("f13Pla1TL",&f13Pla1TL);
  tree->SetBranchAddress("f13Pla1TR",&f13Pla1TR);
  tree->SetBranchAddress("f13Pla2TL",&f13Pla2TL);
  tree->SetBranchAddress("f13Pla2TR",&f13Pla2TR);
	
  tree->SetBranchAddress("beamAoZ",&beamAoZ);
  tree->SetBranchAddress("beamZ",&beamZ);
  
  tree->SetBranchAddress("nebulaNum",&nebulaNum);
  tree->SetBranchAddress("nebulaID",nebulaID);
  tree->SetBranchAddress("nebulaTURaw",nebulaTURaw);
  tree->SetBranchAddress("nebulaTDRaw",nebulaTDRaw);
  tree->SetBranchAddress("nebulaTU",nebulaTU);
  tree->SetBranchAddress("nebulaTU1",nebulaTU1);
  tree->SetBranchAddress("nebulaTU2",nebulaTU2);
  tree->SetBranchAddress("nebulaTD",nebulaTD);
  tree->SetBranchAddress("nebulaTD1",nebulaTD1);
  tree->SetBranchAddress("nebulaTD2",nebulaTD2);
  tree->SetBranchAddress("nebulaTA",nebulaTA);
  tree->SetBranchAddress("nebulaTA1",nebulaTA1);
  tree->SetBranchAddress("nebulaTA2",nebulaTA2);
  tree->SetBranchAddress("nebulaQUPed",nebulaQUPed);
  tree->SetBranchAddress("nebulaQDPed",nebulaQDPed);
  tree->SetBranchAddress("nebulaQPed",nebulaQPed);

  
  TFile *file2 = new TFile(Form("../../root/run%04d.root.NEBslew",runNum),"RECREATE");
  TTree *tree2 = new TTree("NEBslew","NEBslew");
  tree2->Branch("RunNum",&RunNum,"RunNum/I");
  tree2->Branch("EventNum",&EventNum,"EventNum/I");
  tree2->Branch("tofc713",&tofc713,"tofc713/D");
  tree2->Branch("beamAoZ",&beamAoZ,"beamAoZ/D");
  tree2->Branch("beamZ",&beamZ,"beamZ/D");

  tree2->Branch("nebulaNum",&nebulaNum,"nebulaNum/I");
  tree2->Branch("nebulaID",nebulaID,"nebulaID[nebulaNum]/I");
  tree2->Branch("nebulaTURaw",nebulaTURaw,"nebulaTURaw[nebulaNum]/D");
  tree2->Branch("nebulaTDRaw",nebulaTDRaw,"nebulaTDRaw[nebulaNum]/D");
  tree2->Branch("nebulaTU",nebulaTU,"nebulaTU[nebulaNum]/D");
  tree2->Branch("nebulaTU1",nebulaTU1,"nebulaTU1[nebulaNum]/D");
  tree2->Branch("nebulaTU2",nebulaTU2,"nebulaTU2[nebulaNum]/D");
  tree2->Branch("nebulaTD",nebulaTD,"nebulaTD[nebulaNum]/D");
  tree2->Branch("nebulaTD1",nebulaTD1,"nebulaTD1[nebulaNum]/D");
  tree2->Branch("nebulaTD2",nebulaTD2,"nebulaTD2[nebulaNum]/D");
  tree2->Branch("nebulaTA",nebulaTA,"nebulaTA[nebulaNum]/D");
  tree2->Branch("nebulaTA1",nebulaTA1,"nebulaTA1[nebulaNum]/D");
  tree2->Branch("nebulaTA2",nebulaTA2,"nebulaTA2[nebulaNum]/D");
  tree2->Branch("nebulaQUPed",nebulaQUPed,"nebulaQUPed[nebulaNum]/D");
  tree2->Branch("nebulaQDPed",nebulaQDPed,"nebulaQDPed[nebulaNum]/D");
  tree2->Branch("nebulaQPed",nebulaQPed,"nebulaQPed[nebulaNum]/D");

  tree2->Branch("nebulaTOF",nebulaTOF,"nebulaTOF[nebulaNum]/D");
  
  for(Int_t i=0;i<tree->GetEntries();i++)
    {
      if(!(i%1000))
        {
          std::cout<<"Progress rate(%): " << (Double_t)i/(Double_t)tree->GetEntries()*100 << "\r";
          std::cout.flush();
        }

      tree->GetEntry(i);
      if(tree->GetEntryWithIndex(EventNum)<0) continue;

      Double_t t13 = (f13Pla1TL + f13Pla2TL + f13Pla1TR + f13Pla2TR)/4.;

      
      for(Int_t j=0;j<nebulaNum;j++)
	{
	  /*
	    if(nebulaTU[j]<5&&nebulaTD[j]>5)
	    nebulaTU[j] = nebulaTD[j];
	    else if(nebulaTU[j]>5&&nebulaTD[j]<5)
	    nebulaTD[j] = nebulaTU[j];

	    nebulaTA[j] = (nebulaTU[j]+nebulaTD[j])/2.;
	  */	     
	  // nebulaTU[j]  = nebulaTU[j]  - t13;
	  // nebulaTU1[j] = nebulaTU1[j] - t13;
	  // nebulaTU2[j] = nebulaTU2[j] - t13;
	  // nebulaTD[j]  = nebulaTD[j]  - t13;
	  // nebulaTD1[j] = nebulaTD1[j] - t13;
	  // nebulaTD2[j] = nebulaTD2[j] - t13;
	  nebulaTOF[j] = nebulaTA[j] - t13;
	}
      
      tree2->Fill();
    }

  delete nebulaID;
  delete nebulaTURaw;
  delete nebulaTDRaw;
  delete nebulaTU;
  delete nebulaTU1;
  delete nebulaTU2;
  delete nebulaTD;
  delete nebulaTD1;
  delete nebulaTD2;
  delete nebulaTA;
  delete nebulaTA1;
  delete nebulaTA2;
  delete nebulaQUPed;
  delete nebulaQDPed;
  delete nebulaQPed;
  delete nebulaTOF;

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
