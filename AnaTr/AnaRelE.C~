#include "AnaRelE.H"

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCut.h"

void AnaRelE(Int_t runNum)
{
  TFile *file = new TFile(Form("../../root/run%04d.root.FragPID",runNum),"READ");
  TTree *tree;
  file->GetObject("FragPID",tree);
  tree->BuildIndex("EventNum");
  tree->AddFriend("BDC",Form("../../root/run%04d.root.BDC",runNum));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  tree->AddFriend("FDC",Form("../../root/run%04d.root.FDC",runNum));
  tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->AddFriend("Coin",Form("../../root/run%04d.root.Coin",runNum));
  tree->GetFriend("Coin")->BuildIndex("EventNum");
  tree->AddFriend("Neut",Form("../../root/run%04d.root.Neut",runNum));
  tree->GetFriend("Neut")->BuildIndex("EventNum");

  Int_t RunNum, EventNum;
  Double_t t13;
  Double_t bdc1X, bdc1Y, bdc2X, bdc2Y;
  Double_t tgtX,tgtY,tgtZ;
  Double_t fdc1X, fdc1Y, fdc2X, fdc2Y;
  Double_t fdc1A, fdc1B;
  Double_t beamAoZ, beamZ;
  Double_t fragAoZ, fragZ;
  Double_t brho;
  Int_t coinTrigger;
  
  Int_t neutNum;

  Int_t *neutIndex;
  Int_t *neutID;
  Double_t *neutTA;
  Double_t *neutX;
  Double_t *neutY;
  Double_t *neutZ;
  Bool_t *neutNEB;
  Bool_t *neutVETO;
  Double_t *neutFL1;
  Double_t *neutTOF1;
  Double_t *neutBeta1;
  
  neutIndex = new Int_t[200];
  neutID = new Int_t[200];
  neutX = new Double_t[200];
  neutY = new Double_t[200];
  neutZ = new Double_t[200];
  neutTA = new Double_t[200];
  neutNEB = new Bool_t[200];
  neutVETO = new Bool_t[200];
  neutFL1 = new Double_t[200];
  neutTOF1 = new Double_t[200];
  neutBeta1 = new Double_t[200];
  

  tree->SetBranchAddress("RunNum",&RunNum);
  tree->SetBranchAddress("EventNum",&EventNum);

  tree->SetBranchAddress("t13",&t13);
  tree->SetBranchAddress("bdc1X",&bdc1X);
  tree->SetBranchAddress("bdc1Y",&bdc1Y);
  tree->SetBranchAddress("bdc2X",&bdc2X);
  tree->SetBranchAddress("bdc2Y",&bdc2Y);

  tree->SetBranchAddress("tgtX",&tgtX);
  tree->SetBranchAddress("tgtY",&tgtY);
  tree->SetBranchAddress("tgtZ",&tgtZ);
  
  tree->SetBranchAddress("fdc1X",&fdc1X);
  tree->SetBranchAddress("fdc1Y",&fdc1Y);
  tree->SetBranchAddress("fdc2X",&fdc2X);
  tree->SetBranchAddress("fdc2Y",&fdc2Y);
  
  tree->SetBranchAddress("fdc1A",&fdc1A);
  tree->SetBranchAddress("fdc1B",&fdc1B);
  
  tree->SetBranchAddress("beamAoZ",&beamAoZ);
  tree->SetBranchAddress("beamZ",&beamZ);
  
  tree->SetBranchAddress("fragAoZ",&fragAoZ);
  tree->SetBranchAddress("fragZ",&fragZ);
  tree->SetBranchAddress("brho",&brho);

  tree->SetBranchAddress("coinTrigger",&coinTrigger);

  tree->SetBranchAddress("neutNum",&neutNum);
  
  tree->SetBranchAddress("neutID",neutID);
  tree->SetBranchAddress("neutTA",neutTA);
  tree->SetBranchAddress("neutX",neutX);
  tree->SetBranchAddress("neutY",neutY);
  tree->SetBranchAddress("neutZ",neutZ);
  tree->SetBranchAddress("neutNEB",neutNEB);
  tree->SetBranchAddress("neutVETO",neutVETO);
  

  
  Double_t fdc1Z = -2888.82;
  Double_t fdc1a, fdc1b;
  Double_t fdc1ex, fdc1ey, fdc1ez;

  Double_t fragP,fragPx,fragPy,fragPz;
  Double_t fragE,fragKE;

  Double_t neutXX, neutYY, neutZZ, neutTT;

  Double_t neutFL,neutTOF,neutKE,neutP;
  Double_t neutBeta, neutGamma, neutE;
  Double_t neutPx,neutPy,neutPz;
  Bool_t neutIsNEB;
  Double_t relE;
  
  TFile *file2 = new TFile(Form("../../root/run%04d.root.RelE",runNum),"RECREATE");
  TTree *tree2 = new TTree("RelE","RelE");
  tree2->Branch("RunNum",&RunNum,"RunNum/I");
  tree2->Branch("EventNum",&EventNum,"EventNum/I");
  tree2->Branch("beamAoZ",&beamAoZ,"beamAoZ/D");
  tree2->Branch("beamZ",&beamZ,"beamZ/D");
  tree2->Branch("fragAoZ",&fragAoZ,"fragAoZ/D");
  tree2->Branch("fragZ",&fragZ,"fragZ/D");
  tree2->Branch("coinTrigger",&coinTrigger,"coinTrigger/I");
  tree2->Branch("brho",&brho,"brho/D");
  tree2->Branch("fragP",&fragP,"fragP/D");
  tree2->Branch("fragE",&fragE,"fragE/D");
  tree2->Branch("fragKE",&fragKE,"fragKE/D");
  tree2->Branch("neutFL",&neutFL,"neutFL/D");
  tree2->Branch("neutTOF",&neutTOF,"neutTOF/D");
  tree2->Branch("neutBeta",&neutBeta,"neutBeta/D");
  tree2->Branch("neutIsNEB",&neutIsNEB,"neutIsNEB/O");
  tree2->Branch("neutKE",&neutKE,"neutKE/D");
  tree2->Branch("neutP",&neutP,"neutP/D");
  tree2->Branch("relE",&relE,"relE/D");

  tree2->Branch("neutNum",&neutNum,"neutNum/I");
  tree2->Branch("neutX",neutX,"neutX[neutNum]/D");
  tree2->Branch("neutY",neutY,"neutY[neutNum]/D");
  tree2->Branch("neutXX",&neutXX,"neutXX/D");
  tree2->Branch("neutYY",&neutYY,"neutYY/D");
  tree2->Branch("neutZZ",&neutZZ,"neutZZ/D");
  tree2->Branch("neutZ",neutZ,"neutZ[neutNum]/D");
  tree2->Branch("neutTA",neutTA,"neutTA[neutNum]/D");
  tree2->Branch("neutTT",&neutTT,"neutTT/D");
  tree2->Branch("neutNEB",neutNEB,"neutNEB[neutNum]/O");
  tree2->Branch("neutVETO",neutVETO,"neutVETO[neutNum]/O");

  tree2->Branch("neutBeta1",neutBeta1,"neutBeta1[neutNum]/D");
  
  for(Int_t i=0;i<tree->GetEntries();i++)
    {
      if(!(i%1000))
        {
          std::cout<<"Progress rate(%): " << (Double_t)i/(Double_t)tree->GetEntries()*100 << "\r";
          std::cout.flush();
        }

      tree->GetEntry(i);
      if(tree->GetEntryWithIndex(EventNum)<0) continue;


      fdc1a = (fdc1X-tgtX)/(fdc1Z-tgtZ);
      fdc1b = (fdc1Y-tgtY)/(fdc1Z-tgtZ);

      Double_t frac = TMath::Sqrt(fdc1a*fdc1a + fdc1b*fdc1b + 1);
      fdc1ex = fdc1a/frac;
      fdc1ey = fdc1b/frac;
      fdc1ez = 1./frac;
      
      //tentative 32Ne
      Double_t zz = 10;
      Double_t m32ne = 29839.7;
      Double_t amu32ne = 32.0342;
      
      fragP = brho*299.792458*zz;
      fragPx = fragP * fdc1ex;
      fragPy = fragP * fdc1ey;
      fragPz = fragP * fdc1ez;
      fragE = TMath::Sqrt(fragP*fragP + m32ne*m32ne);
      fragKE = (fragE - m32ne)/amu32ne;
      //tentative

      Double_t m_neut = 939.565;
      Double_t amu_neut = 1.00866; 


      for(Int_t j=0;j<neutNum;j++)
	{
	  neutFL1[j] = TMath::Sqrt(TMath::Power(neutX[j]-tgtX,2) + TMath::Power(neutY[j]-tgtY,2) +TMath::Power(neutZ[j],2) );
	  neutTOF1[j] = neutTA[j] - t13;
	  neutBeta1[j] = neutFL1[j]/neutTOF1[j]/299.792458;
	  if(neutBeta1[j]<=0 || TMath::IsNaN(neutBeta1[j]))
	    {
	      neutBeta1[j] = 0;
	      neutTOF1[j] = 0;
	    }
	}

      neutBeta = TMath::MaxElement(neutNum,neutBeta1);

      for(Int_t j=0;j<neutNum;j++)
	{
	  if(neutBeta == neutBeta1[j])
	    {
	      neutFL = neutFL1[j];
	      neutTOF = neutTOF1[j];
	      neutXX = neutX[j];
	      neutYY = neutY[j];
	      neutZZ = neutZ[j];
	      neutTT = neutTA[j];
	      neutIsNEB = neutNEB[j];
	      break;
	    }
	}
      
      
      neutGamma = 1/TMath::Sqrt(1-neutBeta*neutBeta);
      neutE = neutGamma * m_neut;
      neutKE = (neutGamma - 1) * m_neut / amu_neut;
      neutP = neutGamma * m_neut * neutBeta;
      neutPx = neutP *(neutXX-tgtX)/neutFL;
      neutPy = neutP *(neutYY-tgtY)/neutFL;
      neutPz = neutP *(neutZZ)/neutFL;
      

      relE  = TMath::Power(fragE+neutE,2);
      relE -= TMath::Power(fragPx+neutPx,2);
      relE -= TMath::Power(fragPy+neutPy,2);
      relE -= TMath::Power(fragPz+neutPz,2);
      relE = TMath::Sqrt(relE) - m_neut - m32ne;
      
      
      tree2->Fill();
    }

  
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
  
  
  tree2->Write();
  file2->Close();
}
