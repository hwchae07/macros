#include "AnaFragPID.H"

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCut.h"
#include "TEventList.h"

void AnaFragPID(Int_t runNum)
{
  TFile *file = new TFile(Form("../../root/run%04d.root.BeamPID",runNum),"READ");
  TTree *tree;
  file->GetObject("BeamPID",tree);
  tree->BuildIndex("EventNum");
  tree->AddFriend("BDC",Form("../../root/run%04d.root.BDC",runNum));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  tree->AddFriend("FDC",Form("../../root/run%04d.root.FDC",runNum));
  tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->AddFriend("HOD",Form("../../root/run%04d.root.HOD",runNum));
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  tree->AddFriend("SAMURAI",Form("../../root/run%04d.root.SAMURAI",runNum));
  tree->GetFriend("SAMURAI")->BuildIndex("EventNum");
  tree->AddFriend("Coin",Form("../../root/run%04d.root.Coin",runNum));
  tree->GetFriend("Coin")->BuildIndex("EventNum");

  Int_t RunNum,EventNum;
  Double_t f13Pla1TL, f13Pla1TR, f13Pla2TL, f13Pla2TR;
  Double_t betaF7, beamMomU, beamZ, beamAoZ;
  Double_t fdc1X, fdc1Y, fdc1A, fdc1B;
  Double_t fdc2X, fdc2Y, fdc2A, fdc2B;
  Double_t bdc1X, bdc1Y, bdc2X, bdc2Y;
  Int_t hodNum;
  Int_t *hodID;
  Double_t *hodTA;
  Double_t *hodQ;
  Double_t brho;
  Int_t coinTrigger;

  hodID = new Int_t[24];
  hodTA = new Double_t[24];
  hodQ = new Double_t[24];
  
  tree->SetBranchAddress("RunNum",&RunNum);
  tree->SetBranchAddress("EventNum",&EventNum);
  tree->SetBranchAddress("f13Pla1TL",&f13Pla1TL);
  tree->SetBranchAddress("f13Pla1TR",&f13Pla1TR);
  tree->SetBranchAddress("f13Pla2TL",&f13Pla2TL);
  tree->SetBranchAddress("f13Pla2TR",&f13Pla2TR);

  tree->SetBranchAddress("betaF7",&betaF7);
  tree->SetBranchAddress("beamMomU",&beamMomU);
  tree->SetBranchAddress("beamZ",&beamZ);
  tree->SetBranchAddress("beamAoZ",&beamAoZ);

  tree->SetBranchAddress("fdc1X",&fdc1X);
  tree->SetBranchAddress("fdc1A",&fdc1A);
  tree->SetBranchAddress("fdc1Y",&fdc1Y);
  tree->SetBranchAddress("fdc1B",&fdc1B);

  tree->SetBranchAddress("fdc2X",&fdc2X);
  tree->SetBranchAddress("fdc2A",&fdc2A);
  tree->SetBranchAddress("fdc2Y",&fdc2Y);
  tree->SetBranchAddress("fdc2B",&fdc2B);

  tree->SetBranchAddress("bdc1X",&bdc1X);
  tree->SetBranchAddress("bdc1Y",&bdc1Y);
  tree->SetBranchAddress("bdc2X",&bdc2X);
  tree->SetBranchAddress("bdc2Y",&bdc2Y);

  tree->SetBranchAddress("hodNum",&hodNum);
  tree->SetBranchAddress("hodID",hodID);
  tree->SetBranchAddress("hodTA",hodTA);
  tree->SetBranchAddress("hodQ",hodQ);
  
  tree->SetBranchAddress("brho",&brho);
  
  tree->SetBranchAddress("coinTrigger",&coinTrigger);


  Double_t tgtX,tgtY,tgtZ;
  Double_t FL, betaFrag, fragZ, fragAoZ;
  Double_t t13,tof_hod;
  
  TFile *file2 = new TFile(Form("../../root/run%04d.root.FragPID",runNum),"RECREATE");
  TTree *tree2 = new TTree("FragPID","FragPID");
  tree2->Branch("RunNum",&RunNum,"RunNum/I");
  tree2->Branch("EventNum",&EventNum,"EventNum/I");
  tree2->Branch("t13",&t13,"t13/D");
  tree2->Branch("betaF7",&betaF7,"betaF7/D");
  tree2->Branch("beamMomU",&beamMomU,"beamMomU/D");
  tree2->Branch("beamZ",&beamZ,"beamZ/D");
  tree2->Branch("beamAoZ",&beamAoZ,"beamAoZ/D");
  tree2->Branch("tgtX",&tgtX,"tgtX/D");
  tree2->Branch("tgtY",&tgtY,"tgtY/D");
  tree2->Branch("tgtZ",&tgtZ,"tgtZ/D");
  tree2->Branch("brho",&brho,"brho/D");
  tree2->Branch("FL",&FL,"FL/D");
  tree2->Branch("betaFrag",&betaFrag,"betaFrag/D");
  tree2->Branch("fragZ",&fragZ,"fragZ/D");
  tree2->Branch("fragAoZ",&fragAoZ,"fragAoZ/D");
  tree2->Branch("coinTrigger",&coinTrigger,"coinTrigger/I");

  //tree2->Branch("hodxx",&hodxx,"hodxx/D");
  //tree2->Branch("hodzz",&hodzz,"hodzz/D");
  //tree2->Branch("FL1",&FL1,"FL1/D");
  //tree2->Branch("FL2",&FL2,"FL2/D");
  //tree2->Branch("theta",&theta,"theta/D");


  Double_t bdcWidth;
  Double_t bdc2Tgt;
  Double_t bdc1Tgt;




  Double_t fdc1Z, alpha, fdc2Z, fdc2zz, fdc2xx;

  Double_t a1, b1, a2, b2, interZ, interX;

  Double_t hodZ, hodX, hodzz, hodxx, theta;

  Double_t FL1, FL2;

  Double_t dEFacFrag;
  Double_t fragMomUfromHOD;

  
  
  for(Int_t i=0;i<tree->GetEntries();i++)
    {
      if(!(i%1000))
        {
          std::cout<<"Progress rate(%): " << (Double_t)i/(Double_t)tree->GetEntries()*100 << "\r";
          std::cout.flush();
        }

      tree->GetEntry(i);
      if(tree->GetEntryWithIndex(EventNum)<0) continue;

  
      //Target size//
      bdcWidth = 120.;
      bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
      bdc1Tgt = bdc2Tgt + 1000;
      tgtZ = -4978.89+488;

      tgtX = (bdc2X - bdc1X)/1000*bdc1Tgt + bdc1X;
      tgtY = (bdc2Y - bdc1Y)/1000*bdc1Tgt + bdc1Y;
      //Target size//         

      t13 = (f13Pla1TL+f13Pla1TR+f13Pla2TL+f13Pla2TR)/4;
      tof_hod = hodTA[0] - t13;

      if(RunNum>294&&RunNum<=301)
	tof_hod = hodTA[0] - t13 - 2;
      
      //FDC z position//                                                                                                          
      fdc1Z = -2888.82;
      alpha = -59.92/180.*TMath::Pi();
      fdc2Z = 4164.51;    //3726.51+438;
      
      fdc2zz = fdc2Z*TMath::Cos(alpha) - (fdc2X+715.01)*TMath::Sin(alpha);
      fdc2xx = fdc2Z*TMath::Sin(alpha) + (fdc2X+715.01)*TMath::Cos(alpha);
      //FDC z position//                                                                                                          

      //intersection//
      a1 = TMath::Tan(fdc1A);
      b1 = fdc1X - fdc1Z*a1;
      a2 = TMath::Tan(fdc2A + alpha);
      b2 = fdc2xx - fdc2zz*a2;
      interZ = (b2-b1)/(a1-a2);
      interX = (a1*b2-b1*a2)/(a1-a2);
      //intersection//                                                                                                            

      hodZ = 4999.17;

      hodX = fdc2X + 715.01 + (hodZ-fdc2Z)*TMath::Tan(fdc2A);
      hodzz = hodZ*TMath::Cos(alpha) - (hodX)*TMath::Sin(alpha);
      hodxx = hodZ*TMath::Sin(alpha) + (hodX)*TMath::Cos(alpha);
      theta = TMath::Abs(alpha + fdc1A - fdc2A);

      FL1 = sqrt( TMath::Power((interX-tgtX),2)  + TMath::Power((interZ-tgtZ),2) );
      FL2 = sqrt( TMath::Power((interX-hodxx),2) + TMath::Power((interZ-hodzz),2) );
      FL = FL1 + FL2 + 1000 * (brho/2.9) * theta - 2 * 1000 * (brho/2.9) * TMath::Tan(theta/2);

      betaFrag = FL/tof_hod/299.792458;

      dEFacFrag = TMath::Log(15795.98 * betaFrag * betaFrag) - TMath::Log(1-betaFrag*betaFrag) - betaFrag*betaFrag;

      fragZ = TMath::Sqrt(hodQ[0]/TMath::Sqrt( (1+TMath::Power(TMath::Tan(fdc2A),2) + TMath::Power(TMath::Tan(fdc2B),2)) )/dEFacFrag)*betaFrag;
      fragMomUfromHOD = 931.494*betaFrag/TMath::Sqrt(1-betaFrag*betaFrag);
      fragAoZ = brho * 299.792458 / fragMomUfromHOD;


      tree2->Fill();
    }
  tree2->Write();
  file2->Close();

  delete hodID;
  delete hodTA;
  delete hodQ;
  
}
