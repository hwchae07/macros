#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TMath.h"
#include "TH1.h"
#include "TVectorD.h"
#include "TLorentzVector.h"

using namespace TMath;

bool IsFileExist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

TCutG* Load2DGate(const char* filename){

  double* x = new double[50];
  double* y = new double[50];
  int nPoint;
  char* cutName = new char[30];

  ifstream fin(filename);
  fin >> cutName;
  fin >> nPoint;
  if (nPoint == 0)
    return NULL; // No Gate!

  for (int i = 0 ; i < nPoint ; i++)
    fin >> x[i] >> y[i];
  fin.close();

  return (new TCutG(cutName,nPoint,x,y));
}

void Analysis(Int_t runNum, const char* frag) {
  if (IsFileExist(Form("./root/run%04d.root.ADA.%s",runNum,frag))){
    return;}
  
  TFile* file = new TFile(Form("root/run%04d.root.Coin",runNum),"READ");
  TTree *tree = (TTree*)file->Get("Coin");

  tree->BuildIndex("EventNum");
  tree->AddFriend("BeamPla",Form("./root/run%04d.root.BeamPla",runNum));
  tree->AddFriend("BDC",Form("./root/run%04d.root.BDC",runNum));
  //  tree->AddFriend("Neut",Form("./root/run%04d.root.Neut",runNum));
  tree->AddFriend("Brho",Form("./root/run%04d.root.Brho",runNum));
  tree->AddFriend("FDC",Form("./root/run%04d.root.FDC",runNum));
  tree->AddFriend("HOD",Form("./root/run%04d.root.HOD",runNum));
  tree->AddFriend("MINOS",Form("./root/run%04d.root.MINOS",runNum));
  tree->AddFriend("DALI",Form("./root/run%04d.root.DALI",runNum));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  //  tree->GetFriend("Neut")->BuildIndex("EventNum");
  tree->GetFriend("Brho")->BuildIndex("EventNum");
  tree->GetFriend("BeamPla")->BuildIndex("EventNum");
  tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  tree->GetFriend("MINOS")->BuildIndex("EventNum");
  tree->GetFriend("DALI")->BuildIndex("EventNum");

  Int_t EventNum;

  Int_t coinTrigger;

  Double_t f13Pla1TL, f13Pla1TR, f13Pla2TL, f13Pla2TR;
  Double_t f7PlaTL, f7PlaTR, f5PlaTL, f5PlaTR;
  Double_t f13Pla1QL, f13Pla1QR; // Raw
  Double_t f13Pla2QL, f13Pla2QR; // Raw
  Double_t f5X;

  Double_t bdc1X, bdc1Y;
  Double_t bdc2X, bdc2Y;

  Int_t hodNum; // num Pla
  Int_t hodID[24];//[24]; // ID
  Double_t hodQU[24];//[24]; // Cal
  Double_t hodQD[24];//[24]; // Cal
  Double_t hodTA[24];//[24]; // Slw

  Double_t fdc1X, fdc1Y, fdc1A, fdc1B;
  Double_t fdc2X, fdc2Y, fdc2A, fdc2B;

  Double_t fragBrho;

  Int_t neutNum;
  Int_t neulandNum;
  Int_t nebulaNum;
  Int_t neulandMult[9];
  Int_t nebulaMult[6];

  Int_t* neutID;//[200];
  Double_t* neutQU;//[200]; // Cal
  Double_t* neutQD;//[200]; // Cal
  Double_t* neutTU;//[200]; // Cal
  Double_t* neutTD;//[200]; // Cal
  Double_t* neutTA;//[200]; // Slw
  Double_t* neutX;//[200];
  Double_t* neutY;//[200];
  Double_t* neutZ;//[200];
  Bool_t* neutNEB;//[200];
  Bool_t* neutVETO;//[200];
  
  neutID = new Int_t[200];
  neutQU = new Double_t[200];
  neutQD = new Double_t[200];
  neutTU = new Double_t[200];
  neutTD = new Double_t[200];
  neutTA = new Double_t[200];
  neutX = new Double_t[200];
  neutY = new Double_t[200];
  neutZ = new Double_t[200];
  neutNEB = new Bool_t[200];
  neutVETO = new Bool_t[200];

  Int_t minosTrack;
  Double_t minosXv, minosYv, minosZv;
  Double_t minosTrX, minosTrA, minosTrY, minosTrB;

  Int_t daliNum;
  Int_t* daliID;//[144];
  Double_t* daliTheta;//[144]; // Cal
  Double_t* daliE;//[144]; // Cal
  Double_t* daliT;//[144]; // Cal
  daliID = new Int_t[140];
  daliTheta = new Double_t[140];
  daliE = new Double_t[140];
  daliT = new Double_t[140];

  tree->SetBranchAddress("EventNum",&EventNum);
  tree->SetBranchAddress("f5PlaTL",&f5PlaTL);
  tree->SetBranchAddress("f5PlaTR",&f5PlaTR);
  tree->SetBranchAddress("f7PlaTL",&f7PlaTL);
  tree->SetBranchAddress("f7PlaTR",&f7PlaTR);
  tree->SetBranchAddress("f13Pla1TL",&f13Pla1TL);
  tree->SetBranchAddress("f13Pla1TR",&f13Pla1TR);
  tree->SetBranchAddress("f13Pla1QL",&f13Pla1QL);
  tree->SetBranchAddress("f13Pla1QR",&f13Pla1QR);
  tree->SetBranchAddress("f13Pla2TL",&f13Pla2TL);
  tree->SetBranchAddress("f13Pla2TR",&f13Pla2TR);
  tree->SetBranchAddress("f13Pla2QL",&f13Pla2QL);
  tree->SetBranchAddress("f13Pla2QR",&f13Pla2QR);
  tree->SetBranchAddress("f5X",&f5X);

  tree->SetBranchAddress("bdc1X",&bdc1X);
  tree->SetBranchAddress("bdc1Y",&bdc1Y);
  tree->SetBranchAddress("bdc2X",&bdc2X);
  tree->SetBranchAddress("bdc2Y",&bdc2Y);

  tree->SetBranchAddress("hodNum",&hodNum);
  tree->SetBranchAddress("hodID",hodID);
  tree->SetBranchAddress("hodQU",hodQU);
  tree->SetBranchAddress("hodQD",hodQD);
  tree->SetBranchAddress("hodTA",hodTA);

  tree->SetBranchAddress("fdc1X",&fdc1X);
  tree->SetBranchAddress("fdc1A",&fdc1A);
  tree->SetBranchAddress("fdc1Y",&fdc1Y);
  tree->SetBranchAddress("fdc1B",&fdc1B);
  tree->SetBranchAddress("fdc2X",&fdc2X);
  tree->SetBranchAddress("fdc2A",&fdc2A);
  tree->SetBranchAddress("fdc2Y",&fdc2Y);
  tree->SetBranchAddress("fdc2B",&fdc2B);

  tree->SetBranchAddress("brho",&fragBrho);
  /*
  tree->SetBranchAddress("neutNum",&neutNum);
  tree->SetBranchAddress("neulandNum",&neulandNum);
  tree->SetBranchAddress("nebulaNum",&nebulaNum);
  tree->SetBranchAddress("neulandMult",neulandMult);
  tree->SetBranchAddress("nebulaMult",nebulaMult);
  tree->SetBranchAddress("neutID",neutID);
  tree->SetBranchAddress("neutQU",neutQU);
  tree->SetBranchAddress("neutQD",neutQD);
  tree->SetBranchAddress("neutTU",neutTU);
  tree->SetBranchAddress("neutTD",neutTD); 
  tree->SetBranchAddress("neutTA",neutTA);
  tree->SetBranchAddress("neutX",neutX);
  tree->SetBranchAddress("neutY",neutY);
  tree->SetBranchAddress("neutZ",neutZ);
  tree->SetBranchAddress("neutNEB",neutNEB);
  tree->SetBranchAddress("neutVETO",neutVETO);
  */
  tree->SetBranchAddress("minosTrack",&minosTrack);
  tree->SetBranchAddress("minosXv",&minosXv);
  tree->SetBranchAddress("minosYv",&minosYv); 
  tree->SetBranchAddress("minosZv",&minosZv);
  tree->SetBranchAddress("minosTrX",&minosTrX);
  tree->SetBranchAddress("minosTrA",&minosTrA);
  tree->SetBranchAddress("minosTrY",&minosTrY);
  tree->SetBranchAddress("minosTrB",&minosTrB);

  tree->SetBranchAddress("daliNum",&daliNum);
  tree->SetBranchAddress("daliID",daliID);
  tree->SetBranchAddress("daliTheta",daliTheta);
  tree->SetBranchAddress("daliE",daliE);
  tree->SetBranchAddress("daliT",daliT);

  // Const
  const Double_t c = 299.8; // [mm/ns]
  const Double_t amu = 931.494; // [MeV]
  const Double_t neutMass = 939.565; // [MeV/c^2]
  Double_t L7_13 = 36617; // [mm]
  Double_t LTg_Al = 469 + 50; // [mm]
  Double_t LBDC1_BDC2 = 999;
  Double_t LBDC2_MINF = 1110;
  Double_t LMINF_MIN0 = 25;
  Double_t LMIN0_MINR = 125;
  Double_t LBDC = 120;
  // Parameters
  Double_t beamTOFOffset = 235.843; // [ns] for 24O beam
  //  Double_t f5Brho = 7.4498; // [Tm] for 24O beam
  Double_t f5Brho = 7.6172; // [Tm] for 29F beam
  Double_t fragCharge = 0;
  if (frag[2] == 'o') fragCharge = 8;
  else if (frag[2] == 'f') fragCharge = 9;
  Double_t fragMass = 0;
  if (strcmp(frag,"22o") == 0) fragMass = 20498.06444;
  else if (strcmp(frag,"23o") == 0) fragMass = 21434.88726;
  else if (strcmp(frag,"24o") == 0) fragMass = 22370.83871;
  else if (strcmp(frag,"25f") == 0) fragMass = 23294.02403;
  else if (strcmp(frag,"26f") == 0) fragMass = 24232.51711;
  else if (strcmp(frag,"27f") == 0) fragMass = 25170.66621;
  Double_t elossFac = 0;
  if (frag[2] == 'o') elossFac = 4.1684;
  else if (frag[2] == 'f') elossFac = 5.3004;
  Double_t betaCorFac = 0;
  if (strcmp(frag,"22o") == 0) betaCorFac = -0.01432;
  else if (strcmp(frag,"23o") == 0) betaCorFac = -0.01432;
  else if (strcmp(frag,"24o") == 0) betaCorFac = -0.01432;
  else if (strcmp(frag,"25f") == 0) betaCorFac = -0.01459;
  else if (strcmp(frag,"26f") == 0) betaCorFac = -0.0189; 
  else if (strcmp(frag,"27f") == 0) betaCorFac = -0.02232;
  betaCorFac = 0.;
  // Cut
  //  Double_t beamAoZL = 2.95;  Double_t beamAoZH = 3.1;  // run 41
  //  Double_t beamZL = 7;  Double_t beamZH = 9; // run 41
  Double_t beamAoZL = 3.18;  Double_t beamAoZH = 3.25;  // run 210
  Double_t beamZL = 7.6;  Double_t beamZH = 9.7; // run 210
  Double_t tgtR = 28; // [mm]
  Double_t hodQL = 6.8; Double_t hodQH = 8.2;
  if (frag[2] == 'o') { hodQL = 6.8; hodQH = 8.2; }
  else if (frag[2] == 'f') { hodQL = 8.2; hodQH = 9.4; }

  TCutG* fragA = Load2DGate(Form("cut/f%s.cut",frag));
  TFile* saveFile = new TFile(Form("./root/run%04d.root.ADA.%s",runNum,frag),"RECREATE");
  TTree* saveTree = new TTree("Erel","Erel");
  TLorentzVector frag4Mom;
  TLorentzVector neut4Mom;
  Double_t eRel;
  Double_t hopeX, hopeY, hopeZ;
  Bool_t neutNEB2;
  Double_t neutQ, neutT;
  Double_t neutX2, neutY2, neutZ2;
  Double_t daliEMax;
  Double_t daliESum;

  saveTree->Branch("EventNum",&EventNum,"EventNum/I");
  saveTree->Branch("frag4Mom",&frag4Mom);
  saveTree->Branch("neut4Mom",&neut4Mom);
  saveTree->Branch("eRel",&eRel,"eRel/D");
  saveTree->Branch("hopeX",&hopeX,"hopeX/D");
  saveTree->Branch("hopeY",&hopeY,"hopeY/D");
  saveTree->Branch("hopeZ",&hopeZ,"hopeZ/D");
  /*  saveTree->Branch("neutNEB",&neutNEB2,"neutNEB/O");
  saveTree->Branch("neutQ",&neutQ,"neutQ/D");
  saveTree->Branch("neutT",&neutT,"neutT/D");
  saveTree->Branch("neutX",&neutX2,"neutX/D");
  saveTree->Branch("neutY",&neutY2,"neutY/D");
  saveTree->Branch("neutZ",&neutZ2,"neutZ/D");*/
  saveTree->Branch("daliESum",&daliESum,"daliESum/D");
  saveTree->Branch("daliEMax",&daliEMax,"daliEMax/D");

  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  // LOOP
  for (Int_t i = 0 ; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    if (tree->GetEntryWithIndex(EventNum) < 0) continue;
    std::cout << "Event Number : " << EventNum << "\r";

    ////////////////////////////////////////////////////////////////
    // Beam
    Double_t T13 = (f13Pla1TL+f13Pla1TR+f13Pla2TL+f13Pla2TR)/4;
    Double_t Q13 = Sqrt(f13Pla1QL*f13Pla1QR);
    Double_t T7 = (f7PlaTL+f7PlaTR)/2;
    Double_t TOF7_13 = T13 - T7 + beamTOFOffset;

    Double_t beamBeta = L7_13/TOF7_13/c;
    Double_t beamGamma = 1/Sqrt(1-beamBeta*beamBeta);
    Double_t beamBrho = Sqrt(2*amu)*f5Brho/(0.1439)*(1+0.7*f5X/3300);
    Double_t zFactor = Log(2*511*1000*beamBeta*beamBeta/64.7*beamGamma*beamGamma)-beamBeta*beamBeta;

    Double_t AoZ = beamBrho/beamBeta/beamGamma/amu;
    Double_t Z = 1./0.35357*(Sqrt(Q13*beamBeta*beamBeta/zFactor)-4.91653)+8;

    // Cut
    if (AoZ < beamAoZL || AoZ > beamAoZH ||
	Z < beamZL || Z > beamZH)
      continue;
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // BDC
    Double_t tgtFX = bdc2X+(LBDC2_MINF-LBDC/2)/LBDC1_BDC2*(bdc2X-bdc1X);
    Double_t tgtFY = bdc2Y+(LBDC2_MINF-LBDC/2)/LBDC1_BDC2*(bdc2Y-bdc1Y);
    Double_t tgtFR = Sqrt(tgtFX*tgtFX+tgtFY*tgtFY);

    // Cut
    if (tgtFR > tgtR)
      continue;
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // MINOS
    Double_t tgt0X = bdc2X+(LBDC2_MINF+LMINF_MIN0-LBDC/2)/LBDC1_BDC2*(bdc2X-bdc1X);
    Double_t tgt0Y = bdc2Y+(LBDC2_MINF+LMINF_MIN0-LBDC/2)/LBDC1_BDC2*(bdc2Y-bdc1Y);
    Double_t tgt0A = (bdc2X-bdc1X)/LBDC1_BDC2;
    Double_t tgt0B = (bdc2Y-bdc1Y)/LBDC1_BDC2;

    Double_t t = tgt0A*minosTrB - tgt0B*minosTrA;
    Double_t A = (tgt0A - minosTrA)*(tgt0A - minosTrA) + (tgt0B - minosTrB)*(tgt0B - minosTrB) + t*t; 
    Double_t u = (tgt0X - minosTrX)*(minosTrA-tgt0A-minosTrB*t) - (tgt0Y - minosTrY)*(minosTrB - tgt0B - minosTrA*t);
    u /= A;
    Double_t v = (minosTrX - tgt0X)*(tgt0B - minosTrB) + (tgt0Y - minosTrY)*(tgt0A - minosTrA);
    v /= A;

    hopeX = tgt0X + u*tgt0A + 0.5*v*(tgt0B - minosTrB);
    hopeY = tgt0Y + u*tgt0B + 0.5*v*(minosTrA - tgt0A);
    hopeZ = u + 0.5*v*t;
    // MINOS cut
    if (hopeZ < -LMINF_MIN0 || hopeZ > LMIN0_MINR) continue;
    if (hopeX*hopeX+hopeY*hopeY > tgtR*tgtR) continue;
    //    if (minosTrack < 2) continue;
    //    hopeX = minosXv;
    //   hopeY = minosYv;
    //    hopeZ = minosZv;
    ////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////
    // Frag. PID
    Double_t hodQA = Sqrt(hodQU[0]*hodQD[0]);

    // Cut
    if (hodQA < hodQL || hodQA > hodQH)
      continue;
    if (!fragA->IsInside(fdc2X,hodTA[0]))
      continue;
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Frag. mom.
    Double_t fragMomA = Sqrt(2*amu)*fragBrho/(0.1439/fragCharge); // [MeV/c]
    // 
    // Energy loss
    Double_t fragE = Sqrt(fragMomA*fragMomA + fragMass*fragMass);
    fragE += elossFac*(LMIN0_MINR - hopeZ);
    fragMomA = Sqrt(fragE*fragE - fragMass*fragMass);
    //
    // Beta correction
    Double_t fragBeta = fragMomA/fragE;
    fragBeta += betaCorFac;
    Double_t fragGamma = 1/Sqrt(1-fragBeta*fragBeta);
    fragMomA = fragGamma*fragBeta*fragMass;
    //    TVector3 fragMom(fdc1A,fdc1B,1);
    //    TVector3 fragMom(fdc1X-minosXv,fdc1Y-minosYv,1571-minosZv+50);
    TVector3 fragMom(fdc1X-hopeX,fdc1Y-hopeY,1571-hopeZ+50);
    fragMom.SetMag(fragMomA);
    frag4Mom.SetVectM(fragMom,fragMass);
    ////////////////////////////////////////////////////////////////

    /*
    ////////////////////////////////////////////////////////////////
    // Neut. mom.
    // Neutron loop 
    Bool_t nebVETO = false;
    Bool_t neuVETO = false;
    Int_t neutIndex = -1;
    for (Int_t j = 0 ; j < neutNum ; j++){
      nebVETO = nebVETO || (neutVETO[j] && neutNEB[j]);
      neuVETO = neuVETO || (neutVETO[j] && !neutNEB[j]); }
    if (nebulaNum > 0 && neulandNum > 0) continue;
    if (nebVETO || neuVETO) continue;
    //    for (Int_t j = 0 ; j < neutNum ; j++){
    //      if (!neutVETO[neutNum] && 
    //	  ( (!nebVETO && neutNEB[neutNum]) ||
    //	    (!neuVETO && !neutNEB[neutNum]))) neutIndex = neutNum;
    //      break;}
    // choose the fastest one
    Double_t time = 999999;
    for (Int_t j = 0 ; j < neutNum ; j++) {
      if (!neutVETO[j] && 
	  neutTA[j] < time &&
	  //	  !neutNEB[j]) 
	  Sqrt(neutQU[j]*neutQD[j]) > 6) {
	neutIndex = j;
	time = neutTA[j];}}
    if (neutIndex < 0) continue; // Veto cut
    neutNEB2 = neutNEB[neutIndex];
    //    
    Double_t neutTOF = neutTA[neutIndex]  - T13;
    //    cout << neutNum << " " << neutIndex << " " << neutX[neutIndex] << " " << neutY[neutIndex] << " " << neutZ[neutIndex] << " " << minosZv << endl;
    //    neutX[neutIndex] -= minosXv;
    neutX[neutIndex] -= hopeX;
    //    neutY[neutIndex] -= minosYv;
    neutY[neutIndex] -= hopeY;
    neutZ[neutIndex] += LTg_Al;
    //    neutZ[neutIndex] -= minosZv;
    neutZ[neutIndex] -= hopeZ;
    Double_t neutFL = Sqrt(neutX[neutIndex]*neutX[neutIndex]+neutY[neutIndex]*neutY[neutIndex]+neutZ[neutIndex]*neutZ[neutIndex]);
    Double_t neutBeta = neutFL/neutTOF/c;
    Double_t neutGamma = 1/Sqrt(1-neutBeta*neutBeta);
    Double_t neutMomA = neutGamma*neutBeta*neutMass; // [MeV/c]
    TVector3 neutMom(neutX[neutIndex],neutY[neutIndex],neutZ[neutIndex]);
    neutMom.SetMag(neutMomA);
    neut4Mom.SetVectM(neutMom,neutMass);
    neutX2 = neutX[neutIndex];
    neutY2 = neutY[neutIndex];
    neutZ2 = neutZ[neutIndex];
    neutQ = Sqrt(neutQU[neutIndex]*neutQD[neutIndex]);
    neutT = neutTOF;
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // Relative Energy
    TLorentzVector rel4Mom;
    rel4Mom += frag4Mom;
    rel4Mom += neut4Mom;
    
    eRel = rel4Mom.Mag() - fragMass - neutMass; // [MeV]*/
    ////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////
    // DALI
    Int_t daliIndex = -1;
    Double_t maxE = 0;
    daliESum = 0;
    for (Int_t j = 0 ; j < daliNum ; j++){
      if (daliE[j] > maxE) {
	daliIndex = j;
	maxE = daliE[j];}
      daliESum += daliE[j];
    }
    Double_t dopCor = (1 - fragBeta*Cos(daliTheta[daliIndex]))*fragGamma;
    daliEMax = daliE[daliIndex] * dopCor;
    daliESum *= dopCor;
    saveTree->Fill();    
  }
  std::cout << std::endl;
  //  eRelH->Draw();
  saveTree->Write();
  //  saveTree->Draw("eRel>>(50,0,5)","","");
  saveFile->Close();
  delete saveFile;
  
  
  delete neutID;
  delete  neutQU; 
  delete    neutQD; 
  delete    neutTU; 
  delete    neutTD; 
  delete    neutTA; 
  delete    neutX;
  delete    neutY;
  delete    neutZ;
  delete    neutNEB;
  delete    neutVETO;
  file->Close();

  delete file;
  delete fragA;

}
