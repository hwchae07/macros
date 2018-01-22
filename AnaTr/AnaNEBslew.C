#include "AnaNEBslew.H"

#include <iostream>
#include <fstream>
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

  std::ifstream fin2("../../dat/nebula_slew_new.dat");
  Double_t nebula_slew_par0[120]={0,};
  Double_t nebula_slew_par1[120]={0,};
  Double_t nebula_slew_par2[120]={0,};
  Int_t nebula_gamma_entry[120]={0,};
  Double_t dummy;
  for(Int_t i=0;i<120;i++)
    fin2>>dummy>>nebula_slew_par0[i]>>nebula_slew_par1[i]>>nebula_slew_par2[i]>>nebula_gamma_entry[i];
  fin2.close();

  std::ifstream tof("../../dat/nebula_gamma_offset.dat");
  Double_t tof_offset[120]={0,};
  for(Int_t i=0;i<120;i++)
    tof>>tof_offset[i];
  tof.close();

  double slew_nebula_up_par0[120];
  double slew_nebula_up_par1[120];
  double slew_nebula_down_par0[120];
  double slew_nebula_down_par1[120];
  std::ifstream slew_nebula_up("../../dat/slew_2nd_nebula_up_pol1_layer0.dat");
  for(int i=0;i<120;i++)
    slew_nebula_up>>slew_nebula_up_par0[i]>>slew_nebula_up_par1[i];
  slew_nebula_up.close();
  std::ifstream slew_nebula_down("../../dat/slew_2nd_nebula_down_pol1_layer0.dat");
  for(int i=0;i<120;i++)
    slew_nebula_down>>slew_nebula_down_par0[i]>>slew_nebula_down_par1[i];
  slew_nebula_down.close();

  double par_ycal0[144];
  double par_ycal1[144];
  std::ifstream ycal("../../dat/ycal_nebula.dat");
  for(int i = 0 ; i<120; i++)
    ycal >> par_ycal0[i] >> par_ycal1[i];
  ycal.close();

  /*
    std::ifstream test1("../../dat/DTOff_NEBULA.dat");
    for(int i = 0 ; i<144; i++)
    ycal >> par_ycal0[i] >> dummy;
    test1.close();

    std::ifstream test2("../../dat/DTCal_NEBULA.dat");
    for(int i = 0 ; i<144; i++)
    ycal >> dummy >> par_ycal1[i];
    test2.close();
  */

  Int_t RunNum,EventNum;
  Double_t f13Pla1TL,f13Pla1TR,f13Pla2TL,f13Pla2TR;
  Double_t tofc713;
  Double_t beamAoZ, beamZ;

  Int_t nebulaNum;
  Int_t* nebulaID;
  Double_t* nebulaTURaw;
  Double_t* nebulaTDRaw;
  Double_t* nebulaTU;
  Double_t* nebulaTU2;
  Double_t* nebulaTD;
  Double_t* nebulaTD2;
  Double_t* nebulaTA;
  Double_t* nebulaTA2;
  Double_t* nebulaQUPed;
  Double_t* nebulaQDPed;
  Double_t* nebulaQPed;
  Double_t* nebulaQA;

  Double_t* nebulaTOF;

  //U = up, G = gamma, P = position correction, S = slew
  Double_t* nebulaTOFU;
  Double_t* nebulaTOFUS;
  Double_t* nebulaTOFUG;
  Double_t* nebulaTOFUGP;
  Double_t* nebulaTOFUGS;
  Double_t* nebulaTOFUGPS;

  //D = down, G = gamma, P = position correction, S = slew
  Double_t* nebulaTOFD;
  Double_t* nebulaTOFDS;
  Double_t* nebulaTOFDG;
  Double_t* nebulaTOFDGP;
  Double_t* nebulaTOFDGS;
  Double_t* nebulaTOFDGPS;


  Double_t* nebulaDT;
  Double_t* nebulaDTS;

  Double_t* nebulaX;
  Double_t* nebulaY;
  Double_t* nebulaZ;
  Double_t* nebulaY2;

  nebulaID = new Int_t[144];
  nebulaTURaw = new Double_t[144];
  nebulaTDRaw = new Double_t[144];
  nebulaTU = new Double_t[144];
  nebulaTU2 = new Double_t[144];
  nebulaTD = new Double_t[144];
  nebulaTD2 = new Double_t[144];
  nebulaTA = new Double_t[144];
  nebulaTA2 = new Double_t[144];
  nebulaQUPed = new Double_t[144];
  nebulaQDPed = new Double_t[144];
  nebulaQPed = new Double_t[144];
  nebulaQA = new Double_t[144];
  nebulaTOF = new Double_t[144];

  nebulaTOFU    = new Double_t[144];
  nebulaTOFUS   = new Double_t[144];
  nebulaTOFUG   = new Double_t[144];
  nebulaTOFUGP  = new Double_t[144];
  nebulaTOFUGS  = new Double_t[144];
  nebulaTOFUGPS = new Double_t[144];


  nebulaTOFD    = new Double_t[144];
  nebulaTOFDS   = new Double_t[144];
  nebulaTOFDG   = new Double_t[144];
  nebulaTOFDGP  = new Double_t[144];
  nebulaTOFDGS  = new Double_t[144];
  nebulaTOFDGPS = new Double_t[144];


  nebulaDT = new Double_t[144];
  nebulaDTS = new Double_t[144];

  nebulaX = new Double_t[144];
  nebulaY = new Double_t[144];
  nebulaZ = new Double_t[144];
  nebulaY2 = new Double_t[144];

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
  tree->SetBranchAddress("nebulaTU2",nebulaTU2);
  tree->SetBranchAddress("nebulaTD",nebulaTD);
  tree->SetBranchAddress("nebulaTD2",nebulaTD2);
  tree->SetBranchAddress("nebulaTA",nebulaTA);
  tree->SetBranchAddress("nebulaTA2",nebulaTA2);
  tree->SetBranchAddress("nebulaQUPed",nebulaQUPed);
  tree->SetBranchAddress("nebulaQDPed",nebulaQDPed);
  tree->SetBranchAddress("nebulaQPed",nebulaQPed);
  tree->SetBranchAddress("nebulaQA",nebulaQA);

  tree->SetBranchAddress("nebulaDT",nebulaDT);

  tree->SetBranchAddress("nebulaX",nebulaX);
  tree->SetBranchAddress("nebulaY",nebulaY);
  tree->SetBranchAddress("nebulaZ",nebulaZ);

  //tree->SetBranchAddress("nebulaY2",nebulaY2);

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
  tree2->Branch("nebulaTU2",nebulaTU2,"nebulaTU2[nebulaNum]/D");
  tree2->Branch("nebulaTD",nebulaTD,"nebulaTD[nebulaNum]/D");
  tree2->Branch("nebulaTD2",nebulaTD2,"nebulaTD2[nebulaNum]/D");
  tree2->Branch("nebulaTA",nebulaTA,"nebulaTA[nebulaNum]/D");
  tree2->Branch("nebulaTA2",nebulaTA2,"nebulaTA2[nebulaNum]/D");
  tree2->Branch("nebulaQUPed",nebulaQUPed,"nebulaQUPed[nebulaNum]/D");
  tree2->Branch("nebulaQDPed",nebulaQDPed,"nebulaQDPed[nebulaNum]/D");
  tree2->Branch("nebulaQPed",nebulaQPed,"nebulaQPed[nebulaNum]/D");
  tree2->Branch("nebulaQA",nebulaQA,"nebulaQA[nebulaNum]/D");

  tree2->Branch("nebulaTOF",nebulaTOF,"nebulaTOF[nebulaNum]/D");

  tree2->Branch("nebulaTOFU",nebulaTOFU,"nebulaTOFU[nebulaNum]/D");
  tree2->Branch("nebulaTOFUS",nebulaTOFUS,"nebulaTOFUS[nebulaNum]/D");
  tree2->Branch("nebulaTOFUG",nebulaTOFUG,"nebulaTOFUG[nebulaNum]/D");
  tree2->Branch("nebulaTOFUGP",nebulaTOFUGP,"nebulaTOFUGP[nebulaNum]/D");
  tree2->Branch("nebulaTOFUGS",nebulaTOFUGS,"nebulaTOFUGS[nebulaNum]/D");
  tree2->Branch("nebulaTOFUGPS",nebulaTOFUGPS,"nebulaTOFUGPS[nebulaNum]/D");

  tree2->Branch("nebulaTOFD",nebulaTOFD,"nebulaTOFD[nebulaNum]/D");
  tree2->Branch("nebulaTOFDS",nebulaTOFDS,"nebulaTOFDS[nebulaNum]/D");
  tree2->Branch("nebulaTOFDG",nebulaTOFDG,"nebulaTOFDG[nebulaNum]/D");
  tree2->Branch("nebulaTOFDGP",nebulaTOFDGP,"nebulaTOFDGP[nebulaNum]/D");
  tree2->Branch("nebulaTOFDGS",nebulaTOFDGS,"nebulaTOFDGS[nebulaNum]/D");
  tree2->Branch("nebulaTOFDGPS",nebulaTOFDGPS,"nebulaTOFDGPS[nebulaNum]/D");

  tree2->Branch("nebulaDT",nebulaDT,"nebulaDT[nebulaNum]/D");
  tree2->Branch("nebulaDTS",nebulaDTS,"nebulaDTS[nebulaNum]/D");

  tree2->Branch("nebulaX",nebulaX,"nebulaX[nebulaNum]/D");
  tree2->Branch("nebulaY",nebulaY,"nebulaY[nebulaNum]/D");
  tree2->Branch("nebulaZ",nebulaZ,"nebulaZ[nebulaNum]/D");
  tree2->Branch("nebulaY2",nebulaY2,"nebulaY2[nebulaNum]/D");

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
      Double_t v_scint = 150;


      for(Int_t j=0;j<nebulaNum;j++)
	{
	  //v_scint//
	  if(nebulaID[j]<=120)
	    v_scint = 2*par_ycal1[nebulaID[j]-1];
	  else
	    v_scint = 150;
	  //v_scint//

	  Double_t qped = nebulaQPed[j];

	  if(nebulaTU[j] == -9999 || nebulaTD[j] == -9999)
	    {
	      nebulaTOFUS[j] = TMath::QuietNaN();
	      nebulaTOFDS[j] = TMath::QuietNaN();
	      nebulaY2[j] = TMath::QuietNaN();
	      nebulaY[j] = TMath::QuietNaN();
	    }
	  else
	    {
	      nebulaTOFUS[j]
		= nebulaTU[j] - t13
		- slew_nebula_up_par0[nebulaID[j]-1]
		- slew_nebula_up_par1[nebulaID[j]-1]*TMath::Power(qped,-0.5);

	      nebulaTOFDS[j]
		= nebulaTD[j] - t13
		- slew_nebula_down_par0[nebulaID[j]-1]
		- slew_nebula_down_par1[nebulaID[j]-1] * TMath::Power(qped,-0.5);

	      if(nebulaID[j]<=120)
		{
		  nebulaY2[j]
		    = par_ycal0[nebulaID[j]-1]
		    + par_ycal1[nebulaID[j]-1] * (nebulaTOFDS[j] - nebulaTOFUS[j]);

		  nebulaY[j]
		    = par_ycal0[nebulaID[j]-1]
		    + par_ycal1[nebulaID[j]-1] * (nebulaTOFDS[j] - nebulaTOFUS[j]);
		}

	    }



	  if(nebulaTU[j] == -9999)
	    {
	      nebulaTOFU[j] = TMath::QuietNaN();
	      nebulaTOFUG[j] = TMath::QuietNaN();
	      nebulaTOFUGS[j] = TMath::QuietNaN();
	      nebulaTOFUGP[j] = TMath::QuietNaN();
	      nebulaTOFUGPS[j] = TMath::QuietNaN();
	    }
	  else
	    {
	      nebulaTOFU[j] = nebulaTU[j] - t13;
	      nebulaTOFU[j] = nebulaTOFU[j] - tof_offset[nebulaID[j]-1];

	      double distNEB = TMath::Sqrt(nebulaX[j]*nebulaX[j] + nebulaY[j]*nebulaY[j] + nebulaZ[j]*nebulaZ[j]);
	      nebulaTOFUG[j]  = nebulaTOFU[j] - distNEB/299.792458;

	      nebulaTOFUGS[j]
		= nebulaTOFUG[j]
		- slew_nebula_up_par0[nebulaID[j]-1]
		- slew_nebula_up_par1[nebulaID[j]-1]*TMath::Power(qped,-0.5);

	      nebulaTOFUGP[j] = nebulaTOFUG[j] + nebulaY[j]/v_scint;

	      nebulaTOFUGPS[j]
		= nebulaTOFUGP[j]
		- slew_nebula_up_par0[nebulaID[j]-1]
		- slew_nebula_up_par1[nebulaID[j]-1]*TMath::Power(qped,-0.5);

	    }

	  if(nebulaTD[j] == -9999)
	    {
	      nebulaTOFD[j] = TMath::QuietNaN();
	      nebulaTOFDG[j] = TMath::QuietNaN();
	      nebulaTOFDGS[j] = TMath::QuietNaN();
	      nebulaTOFDGP[j] = TMath::QuietNaN();
	      nebulaTOFDGPS[j] = TMath::QuietNaN();
	    }
	  else
	    {
	      nebulaTOFD[j] = nebulaTD[j] - t13;
	      nebulaTOFD[j] = nebulaTOFD[j] - tof_offset[nebulaID[j]-1];

	      double distNEB = TMath::Sqrt(nebulaX[j]*nebulaX[j] + nebulaY[j]*nebulaY[j] + nebulaZ[j]*nebulaZ[j]);
	      nebulaTOFDG[j] = nebulaTOFD[j] - distNEB/299.792458;

	      nebulaTOFDGS[j]
		= nebulaTOFDG[j]
		- slew_nebula_down_par0[nebulaID[j]-1]
		- slew_nebula_down_par1[nebulaID[j]-1] * TMath::Power(qped,-0.5);

	      nebulaTOFDGP[j] = nebulaTOFDG[j] - nebulaY[j]/v_scint;

	      nebulaTOFDGPS[j]
		= nebulaTOFDGP[j]
		- slew_nebula_down_par0[nebulaID[j]-1]
		- slew_nebula_down_par1[nebulaID[j]-1] * TMath::Power(qped,-0.5);
	    }

	  if(nebulaTU[j] == -9999 || nebulaTD[j] == -9999)
	    {
	      nebulaTOF[j] = TMath::QuietNaN();
	      nebulaDTS[j] = TMath::QuietNaN();
	    }
	  else
	    {
	      nebulaTOF[j] = nebulaTA[j] - t13;
	      nebulaTOF[j] = nebulaTOF[j] - tof_offset[nebulaID[j]-1];
	      double distNEB = TMath::Sqrt(nebulaX[j]*nebulaX[j] + nebulaY[j]*nebulaY[j] + nebulaZ[j]*nebulaZ[j]);
	      nebulaTOF[j] = nebulaTA[j] - distNEB/299.792458;

	      nebulaDTS[j] = nebulaTOFUGS[j] - nebulaTOFDGS[j];
	    }

	}

      tree2->Fill();
    }

  delete nebulaID;
  delete nebulaTURaw;
  delete nebulaTDRaw;
  delete nebulaTU;
  delete nebulaTU2;
  delete nebulaTD;
  delete nebulaTD2;
  delete nebulaTA;
  delete nebulaTA2;
  delete nebulaQUPed;
  delete nebulaQDPed;
  delete nebulaQPed;

  delete nebulaQA;

  delete nebulaTOF;

  delete nebulaTOFU;
  delete nebulaTOFUS;
  delete nebulaTOFUG;
  delete nebulaTOFUGP;
  delete nebulaTOFUGS;
  delete nebulaTOFUGPS;

  delete nebulaTOFD;
  delete nebulaTOFDS;
  delete nebulaTOFDG;
  delete nebulaTOFDGP;
  delete nebulaTOFDGS;
  delete nebulaTOFDGPS;

  delete nebulaDT;
  delete nebulaDTS;

  delete nebulaX;
  delete nebulaY;
  delete nebulaZ;
  delete nebulaY2;

  tree2->Print();
  tree2->Write();
  file2->Close();

}
