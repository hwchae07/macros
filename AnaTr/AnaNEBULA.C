#include <iostream>
#include <fstream>

#include "AnaNEBULA.H"

#include "TTree.h"
#include "TMath.h"

#include "loadPar.H"

AnaNEBULA::AnaNEBULA():
  AnaModule("NEBULA"),
  fNEBULAParameters(NULL),
  fCalibNEBULA(NULL) {
  nebulaNum = 0;
  nebulaID = new Int_t[144];
  nebulaQU = new Double_t[144];
  nebulaQD = new Double_t[144];
  nebulaTU = new Double_t[144];
  nebulaTD = new Double_t[144];
  nebulaTA = new Double_t[144];
  nebulaX = new Double_t[144];
  nebulaY = new Double_t[144];
  //nebulaY2 = new Double_t[144];
  nebulaZ = new Double_t[144];
  nebulaTURaw = new Double_t[144];
  nebulaTDRaw = new Double_t[144];

  nebulaTU2 = new Double_t[144];
  nebulaTD2 = new Double_t[144];
  nebulaTA2 = new Double_t[144];
  nebulaDT = new Double_t[144];

  nebulaQURaw = new Double_t[144];
  nebulaQDRaw = new Double_t[144];
  nebulaQUPed = new Double_t[144];
  nebulaQDPed = new Double_t[144];
  nebulaQPed = new Double_t[144];

  nebulaQA = new Double_t[144];

  //tofOffset = new Double_t[144];
  //par_ycal0 = new double[120];
  //par_ycal1 = new double[120];
  TAveOff = new double[144];

  tdcAlignU = new double[120];
  tdcAlignD = new double[120];
}

AnaNEBULA::~AnaNEBULA(){
  DeleteAll();
  delete nebulaID;
  delete nebulaQU;
  delete nebulaQD;
  delete nebulaTU;
  delete nebulaTD;
  delete nebulaTA;
  delete nebulaX;
  delete nebulaY;
  //delete nebulaY2;
  delete nebulaZ;
  //delete tofOffset;
  delete nebulaTURaw;
  delete nebulaTDRaw;

  delete nebulaTU2;
  delete nebulaTD2;
  delete nebulaTA2;
  delete nebulaDT;

  delete nebulaQURaw;
  delete nebulaQDRaw;
  delete nebulaQUPed;
  delete nebulaQDPed;

  delete nebulaQPed;
  //delete nebulaT;

  delete nebulaQA;

  delete par_ycal0;
  delete par_ycal1;
  delete TAveOff;

  delete tdcAlignU;
  delete tdcAlignD;
}

void AnaNEBULA::InitParameter(){
  fNEBULAParameters = TArtSAMURAIParameters::Instance();
  fNEBULAParameters->LoadParameter((char*)"db/NEBULA.xml");
  /*
    std::ifstream tof("data/nebula_tof_offset.dat");
    Double_t temp;
    for (Int_t i = 0 ; i < 144 ; i++)
    tof >> tofOffset[i] >> temp;

    tof.close();
  */

  // tdc calibration using 6th pol //
  int order = 6;
  parTU = loadParU(order);
  parTD = loadParU(order);
  // tdc calibration using 6th pol //

  /*
    std::ifstream ycal("../../dat/ycal_nebula.dat");
    for(int i = 0 ; i<120; i++)
    ycal >> par_ycal0[i] >> par_ycal1[i];
    ycal.close();
  */

  std::ifstream taveoff("../../dat/NEBULA_TAveOff.dat");
  double dummy;
  for(int i = 0 ; i<144; i++)
    taveoff>> dummy >>TAveOff[i];
  taveoff.close();

  std::ifstream tdcalignu("../../dat/tdc_align_nebula_up.dat");
  for(int i = 0 ; i<120; i++)
    tdcalignu>>tdcAlignU[i];
  tdcalignu.close();

  std::ifstream tdcalignd("../../dat/tdc_align_nebula_down.dat");
  for(int i=0 ;i<120;i++)
    tdcalignd>>tdcAlignD[i];
  tdcalignd.close();


  parLoaded = true;}

void AnaNEBULA::InitDetector(){
  if (!parLoaded) return;
  fCalibNEBULA = new TArtCalibNEBULA;
  detLoaded = true;}

void AnaNEBULA::Analysis(){
  if (!detLoaded) return;
  anaFlag = true;

  fCalibNEBULA->ClearData();
  fCalibNEBULA->ReconstructData();

  nebulaNum = 0;
  for (Int_t i = 0 ; i < 6 ; i++)
    nebulaMult[i] = 0;

  std::ifstream fin("../../dat/nebula_gamma_offset.dat");
  Double_t nebula_offset[120]={0,};
  for(Int_t i=0;i<120;i++)
    fin>>nebula_offset[i];
  fin.close();


  for (Int_t i = 0 ; i < fCalibNEBULA->GetNumNEBULAPla() ; i++)
    {
      TArtNEBULAPla* pla = fCalibNEBULA->GetNEBULAPla(i);
      if (pla->GetHit()){
	nebulaID[nebulaNum] = pla->GetID();
	if (nebulaID[nebulaNum] < 121)
	  nebulaMult[(nebulaID[nebulaNum]-1)/30]++;
	else if (nebulaID[nebulaNum] < 145)
	  nebulaMult[(nebulaID[nebulaNum]-121)/12+4]++;
	else
	  continue;
	nebulaQU[nebulaNum] = pla->GetQUCal();
	nebulaQD[nebulaNum] = pla->GetQDCal();
	nebulaQA[nebulaNum] = TMath::Sqrt(nebulaQU[nebulaNum]*nebulaQD[nebulaNum]);

	double TUCal = pla->GetTURaw();
	double TDCal = pla->GetTDRaw();

	nebulaTU[nebulaNum] = 0;
	nebulaTD[nebulaNum] = 0;

	int order = 6;

	for(int j=0;j<=order;j++)
	  {
	    nebulaTU[nebulaNum] += parTU[nebulaID[nebulaNum]-1][j]*TMath::Power(TUCal,j);
	    nebulaTD[nebulaNum] += parTD[nebulaID[nebulaNum]-1][j]*TMath::Power(TDCal,j);
	  }

	nebulaTU[nebulaNum] += tdcAlignU[nebulaID[nebulaNum]-1];
	nebulaTD[nebulaNum] += tdcAlignD[nebulaID[nebulaNum]-1];

	nebulaTU2[nebulaNum] = pla->GetTUCal();
	nebulaTD2[nebulaNum] = pla->GetTDCal();
	//nebulaTU2[nebulaNum] = pla->GetTUSlw();
	//nebulaTD2[nebulaNum] = pla->GetTDSlw();

	nebulaDT[nebulaNum] = nebulaTU[nebulaNum] - nebulaTD[nebulaNum];

	//nebulaY2[nebulaNum] = par_ycal0[nebulaID[nebulaNum]-1] + par_ycal1[nebulaID[nebulaNum]-1] * (nebulaTD[nebulaNum] - nebulaTU[nebulaNum]);


	if(TUCal<0)
	  {
	    nebulaTU[nebulaNum] = TMath::QuietNaN();
	    nebulaTU2[nebulaNum] = TMath::QuietNaN();
	    nebulaDT[nebulaNum] = TMath::QuietNaN();
	    //nebulaY2[nebulaNum] = TMath::QuietNaN();
	  }
	if(TDCal<0)
	  {
	    nebulaTD[nebulaNum] = TMath::QuietNaN();
	    nebulaTD2[nebulaNum] = TMath::QuietNaN();
	    nebulaDT[nebulaNum] = TMath::QuietNaN();
	    //nebulaY2[nebulaNum] = TMath::QuietNaN();
	  }



	//nebulaTA[nebulaNum] = pla->GetTAveSlw();
	nebulaTA[nebulaNum] = (nebulaTU[nebulaNum] + nebulaTD[nebulaNum])/2.;// + TAveOff[nebulaID[nebulaNum]-1];
	nebulaTA2[nebulaNum] = (nebulaTU2[nebulaNum] + nebulaTD2[nebulaNum])/2.;

	/*
	  nebulaTU[nebulaNum] = pla->GetTUSlw() - nebula_offset[nebulaID[nebulaNum]-1];
	  nebulaTD[nebulaNum] = pla->GetTDSlw() - nebula_offset[nebulaID[nebulaNum]-1];
	  nebulaTA[nebulaNum] = pla->GetTAveSlw() - nebula_offset[nebulaID[nebulaNum]-1];
	*/
	//      std::cout<<nebula_offset[nebulaID[nebulaNum]]<<std::endl;
	nebulaX[nebulaNum] = pla->GetPos(0);

	nebulaY[nebulaNum] = pla->GetPos(1);
	nebulaZ[nebulaNum] = pla->GetPos(2);
	nebulaTURaw[nebulaNum] = pla->GetTURaw();
	nebulaTDRaw[nebulaNum] = pla->GetTDRaw();
	nebulaQURaw[nebulaNum] = pla->GetQURaw();
	nebulaQDRaw[nebulaNum] = pla->GetQDRaw();
	nebulaQUPed[nebulaNum] = pla->GetQUPed();
	nebulaQDPed[nebulaNum] = pla->GetQDPed();

	nebulaQPed[nebulaNum] = TMath::Sqrt(nebulaQUPed[nebulaNum]*nebulaQDPed[nebulaNum]);

	/*
	  Double_t qped = nebulaQPed[nebulaNum];
	  nebulaTUSlew[nebulaNum]
	  = nebulaTU[nebulaNum]
	  - nebula_slew_par1[nebulaID[nebulaNum]-1]*TMath::Power(qped,-1./4.)
	  - nebula_slew_par2[nebulaID[nebulaNum]-1]*TMath::Power(qped,-1./2.);

	  nebulaTDSlew[nebulaNum]
	  = nebulaTD[nebulaNum]
	  - nebula_slew_par1[nebulaID[nebulaNum]-1]*TMath::Power(qped,-1./4.)
	  - nebula_slew_par2[nebulaID[nebulaNum]-1]*TMath::Power(qped,-1./2.);

	  nebulaTASlew[nebulaNum] = (nebulaTUSlew[nebulaNum] + nebulaTDSlew[nebulaNum])/2.;
	*/

	nebulaNum++;
      }
    }
  if (nebulaNum < 1)
    anaFlag = false;

  //delete parTU;
  //delete parTD;
}

void AnaNEBULA::DeleteAll(){
  if (fCalibNEBULA)    delete fCalibNEBULA;
  if (fNEBULAParameters) delete fNEBULAParameters;
  if (parTU) delete parTU;
  if (parTD) delete parTD;
}

void AnaNEBULA::SetTree(){
  AnaModule::SetTree();
  tree->Branch("nebulaNum",&nebulaNum,"nebulaNum/I");
  tree->Branch("nebulaMult",nebulaMult,"nebulaMult[6]/I");
  tree->Branch("nebulaID",nebulaID,"nebulaID[nebulaNum]/I");

  tree->Branch("nebulaQU",nebulaQU,"nebulaQU[nebulaNum]/D");
  tree->Branch("nebulaQD",nebulaQD,"nebulaQD[nebulaNum]/D");
  tree->Branch("nebulaTU",nebulaTU,"nebulaTU[nebulaNum]/D");
  tree->Branch("nebulaTD",nebulaTD,"nebulaTD[nebulaNum]/D");
  tree->Branch("nebulaTA",nebulaTA,"nebulaTA[nebulaNum]/D");
  tree->Branch("nebulaTU2",nebulaTU2,"nebulaTU2[nebulaNum]/D");
  tree->Branch("nebulaTD2",nebulaTD2,"nebulaTD2[nebulaNum]/D");
  tree->Branch("nebulaTA2",nebulaTA2,"nebulaTA2[nebulaNum]/D");
  tree->Branch("nebulaX",nebulaX,"nebulaX[nebulaNum]/D");
  tree->Branch("nebulaY",nebulaY,"nebulaY[nebulaNum]/D");
  //tree->Branch("nebulaY2",nebulaY2,"nebulaY2[nebulaNum]/D");
  tree->Branch("nebulaZ",nebulaZ,"nebulaZ[nebulaNum]/D");
  tree->Branch("nebulaTURaw",nebulaTURaw,"nebulaTURaw[nebulaNum]/D");
  tree->Branch("nebulaTDRaw",nebulaTDRaw,"nebulaTDRaw[nebulaNum]/D");

  tree->Branch("nebulaDT",nebulaDT,"nebulaDT[nebulaNum]/D");

  tree->Branch("nebulaQURaw",nebulaQURaw,"nebulaQURaw[nebulaNum]/D");
  tree->Branch("nebulaQDRaw",nebulaQDRaw,"nebulaQDRaw[nebulaNum]/D");

  tree->Branch("nebulaQUPed",nebulaQUPed,"nebulaQUPed[nebulaNum]/D");
  tree->Branch("nebulaQDPed",nebulaQDPed,"nebulaQDPed[nebulaNum]/D");
  tree->Branch("nebulaQPed",nebulaQPed,"nebulaQPed[nebulaNum]/D");

  tree->Branch("nebulaQA",nebulaQA,"nebulaQA[nebulaNum]/D");



  //tree->SetAlias("nebulaQA","sqrt(nebulaQU*nebulaQD)");
}
