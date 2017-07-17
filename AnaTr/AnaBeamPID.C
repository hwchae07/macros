#include "AnaBeamPID.H"

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCut.h"
#include "TEventList.h"

void AnaBeamPID(Int_t runNum)
{

  TFile *plaTCut = new TFile("../../cut/plaTCut.root","READ");
  TCut *f3pla   = (TCut*)gDirectory->Get("f3plat");
  TCut *f5pla   = (TCut*)gDirectory->Get("f5plat");
  TCut *f7pla   = (TCut*)gDirectory->Get("f7plat");
  TCut *f13pla1 = (TCut*)gDirectory->Get("f13plat1");
  TCut *f13pla2 = (TCut*)gDirectory->Get("f13plat2");

  plaTCut->Close();  
  
  TFile *file = new TFile(Form("../../root/run%04d.root.BeamPla",runNum),"READ");
  TTree *tree;
  file->GetObject("BeamPla",tree);
  tree->BuildIndex("RunNum","EventNum");

  Int_t RunNum,EventNum;
  Double_t f7PlaQL, f7PlaQR;
  Double_t f3PlaTL, f3PlaTR; // Cal
  Double_t f5PlaTL, f5PlaTR; // Cal
  Double_t f7PlaTL, f7PlaTR; // Cal
  Double_t f13Pla1TL, f13Pla1TR; // Cal
  Double_t f13Pla2TL, f13Pla2TR; // Cal
  Double_t f5X;
  Double_t icbEloss;
  
  tree->SetBranchAddress("RunNum",&RunNum);
  tree->SetBranchAddress("EventNum",&EventNum);
  tree->SetBranchAddress("f3PlaTL",&f3PlaTL);
  tree->SetBranchAddress("f3PlaTR",&f3PlaTR);
  tree->SetBranchAddress("f5PlaTL",&f5PlaTL);
  tree->SetBranchAddress("f5PlaTR",&f5PlaTR);
  tree->SetBranchAddress("f7PlaTL",&f7PlaTL);
  tree->SetBranchAddress("f7PlaTR",&f7PlaTR);
  tree->SetBranchAddress("f7PlaQL",&f7PlaQL);
  tree->SetBranchAddress("f7PlaQR",&f7PlaQR);
  tree->SetBranchAddress("f13Pla1TL",&f13Pla1TL);
  tree->SetBranchAddress("f13Pla1TR",&f13Pla1TR);
  tree->SetBranchAddress("f13Pla2TL",&f13Pla2TL);
  tree->SetBranchAddress("f13Pla2TR",&f13Pla2TR);
  tree->SetBranchAddress("f5X",&f5X);
  tree->SetBranchAddress("icbEloss",&icbEloss);

  tree->SetAlias("T3","(f3PlaTL+f3PlaTR)/2");
  tree->SetAlias("T7","(f7PlaTL+f7PlaTR)/2");
  tree->SetAlias("T13","(f13Pla1TL+f13Pla1TR+f13Pla2TL+f13Pla2TR)/4");
  tree->SetAlias("TOF3_13","T13-T3");
  tree->SetAlias("TOF7_13","T13-T7");
  tree->SetAlias("TOF3_7","T7-T3");

  Double_t tof713_offset = 213.73055 + 352.345;   //nsec
  Double_t fl7_13 = (3957.1684 + 3957.2192) / 2.; //cm
  Double_t brho0 = 7.8412; //Tm
  Double_t dispersion = 3280;

  Double_t tofc713;
  Double_t Q7;
  Double_t betaF7, betaF5, betaF13;
  Double_t dEFac_beam,beamMomU;
  Double_t beamZ;
  Double_t beamAoZ;

  TFile *file2 = new TFile(Form("../../root/run%04d.root.BeamPID",runNum),"RECREATE");
  TTree *tree2 = new TTree("BeamPID","BeamPID");
  tree2->Branch("RunNum",&RunNum,"RunNum/I");
  tree2->Branch("EventNum",&EventNum,"EventNum/I");
  tree2->Branch("betaF7",&betaF7,"betaF7/D");
  tree2->Branch("brho0",&brho0,"brho0/D");
  tree2->Branch("f5X",&f5X,"f5X/D");
  tree2->Branch("f7PlaTL",&f7PlaTL,"f7PlaTL/D");
  tree2->Branch("f7PlaTR",&f7PlaTR,"f7PlaTR/D");
  tree2->Branch("f13Pla1TL",&f13Pla1TL,"f13Pla1TL/D");
  tree2->Branch("f13Pla1TR",&f13Pla1TR,"f13Pla1TR/D");
  tree2->Branch("f13Pla2TL",&f13Pla2TL,"f13Pla2TL/D");
  tree2->Branch("f13Pla2TR",&f13Pla2TR,"f13Pla2TR/D");
  tree2->Branch("tofc713",&tofc713,"tofc713/D");
  tree2->Branch("Q7",&Q7,"Q7/D");
  tree2->Branch("beamMomU",&beamMomU,"beamMomU/D");
  tree2->Branch("beamZ",&beamZ,"beamZ/D");
  tree2->Branch("beamAoZ",&beamAoZ,"beamAoZ/D");

  tree2->SetAlias("T3","(f3PlaTL+f3PlaTR)/2");
  tree2->SetAlias("T7","(f7PlaTL+f7PlaTR)/2");
  tree2->SetAlias("T13","(f13Pla1TL+f13Pla1TR+f13Pla2TL+f13Pla2TR)/4");
  tree2->SetAlias("TOF3_13","T13-T3");
  tree2->SetAlias("TOF7_13","T13-T7");
  tree2->SetAlias("TOF3_7","T7-T3");

  
  //plastic TDC cut// 
  tree->Draw(">>evtlist",f3pla[0]&&f5pla[0]&&f7pla[0]&&f13pla1[0]&&f13pla2[0]);
  TEventList *evtlist = (TEventList*)gDirectory->Get("evtlist");
  tree->SetEventList(evtlist);
  //plastic TDC cut// 
  
  for(Int_t i=0 ; i<tree->GetEntries() ; i++)
    {
      if(!(i%1000))
	{
	  std::cout<<"Progress rate(%): " << (Double_t)i/(Double_t)tree->GetEntries()*100 << "\r";
          std::cout.flush();
	}
      
      tree->GetEntry(i);
      if(tree->GetEntryWithIndex(RunNum,EventNum)<0) continue;

      Q7 = TMath::Sqrt(f7PlaQL*f7PlaQR);
      
      Double_t t7  = (f7PlaTL + f7PlaTR) / 2.;
      Double_t t13 = (f13Pla1TL + f13Pla1TR + f13Pla2TL + f13Pla2TR) / 4.;
      tofc713 = (t13 - t7) + tof713_offset;
      betaF7  = fl7_13 / tofc713 / 29.9792458;
      betaF5  = betaF7 * 1.005;
      betaF13 = betaF7 * 0.9;

      
      dEFac_beam = TMath::Log(4866 * betaF13*betaF13) - TMath::Log(1-betaF13*betaF13) - betaF13*betaF13;
      beamZ      = 1+10*TMath::Sqrt(icbEloss/dEFac_beam)*betaF13;
      beamMomU   = 931.494*betaF5/TMath::Sqrt(1-betaF5*betaF5);
      beamAoZ    = brho0 * (1+f5X/dispersion) * 299.792458 / beamMomU;

      if( evtlist->GetIndex(i) != -1)
	tree2->Fill();
      
    }

  tree2->Write();
  file2->Close();

}
