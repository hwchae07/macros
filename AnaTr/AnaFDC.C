#include "AnaFDC.H"

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "TArtDCTrack.hh"

AnaFDC::AnaFDC():
  AnaModule("FDC"),
  fFDCParameters(NULL),
  fCalibFDC1Hit(NULL),
  fCalibFDC1Track(NULL),
  fCalibFDC2Hit(NULL),
  fCalibFDC2Track(NULL) {
  ;}

AnaFDC::~AnaFDC(){    
  DeleteAll();}

void AnaFDC::InitParameter(){
  fFDCParameters = TArtSAMURAIParameters::Instance();
  fFDCParameters->LoadParameter((char*)"db/SAMURAIFDC1.xml");
  fFDCParameters->LoadParameter((char*)"db/SAMURAIFDC2.xml");
  parLoaded = true;}

void AnaFDC::InitDetector(){
  if (!parLoaded) return;
  fCalibFDC1Hit = new TArtCalibFDC1Hit;
  fCalibFDC1Track = new TArtCalibFDC1Track;
  fCalibFDC2Hit = new TArtCalibFDC2Hit;
  fCalibFDC2Track = new TArtCalibFDC2Track;

  LoadTDCDistribution((char*)"db/dc/s027_run0167.root");
  detLoaded = true;}

void AnaFDC::Analysis(){
  anaFlag = true;
  if (!detLoaded) return;

  fCalibFDC1Hit->ClearData();
  fCalibFDC1Track->ClearData();
  fCalibFDC2Hit->ClearData();
  fCalibFDC2Track->ClearData();

  fCalibFDC1Hit->ReconstructData();
  fCalibFDC1Track->ReconstructData();
  fCalibFDC2Hit->ReconstructData();
  fCalibFDC2Track->ReconstructData();

  if (fCalibFDC1Track->GetNumDCTrack() < 1) {
    anaFlag = false; return; }
  if (fCalibFDC2Track->GetNumDCTrack() < 1) {
    anaFlag = false; return; }

  if (fCalibFDC1Track->GetNumDCTrack()) {
    TArtDCTrack *trk = fCalibFDC1Track->GetDCTrack(0);
    Double_t chi2 = trk->GetChi2();
    Int_t ndf = trk->GetNDF();
    Double_t posx = trk->GetPosition(0);
    Double_t posy = trk->GetPosition(1);
    Double_t angx = TMath::ATan(trk->GetAngle(0));
    Double_t angy = TMath::ATan(trk->GetAngle(1));
    fdc1X = posx; fdc1A = angx;
    fdc1Y = posy; fdc1B = angy;
    fdc1C = chi2;

    Double_t frac = TMath::Sqrt(1 + trk->GetAngle(0)*trk->GetAngle(0) + trk->GetAngle(1)*trk->GetAngle(1));
    fdc1ex = trk->GetAngle(0)/frac;
    fdc1ey = trk->GetAngle(1)/frac;
    fdc1ez = 1./frac;

  }

  if (fCalibFDC2Track->GetNumDCTrack()) {
    TArtDCTrack *trk = fCalibFDC2Track->GetDCTrack(0);
    Double_t chi2 = trk->GetChi2();
    Int_t ndf = trk->GetNDF();
    Double_t posx = trk->GetPosition(0);
    Double_t posy = trk->GetPosition(1);
    Double_t angx = TMath::ATan(trk->GetAngle(0));
    Double_t angy = TMath::ATan(trk->GetAngle(1));
    fdc2X = posx; fdc2A = angx;
    fdc2Y = posy; fdc2B = angy;
    fdc2C = chi2;


  }
}

void AnaFDC::DeleteAll(){
  if (fCalibFDC1Hit)    delete fCalibFDC1Hit;
  if (fCalibFDC1Track)    delete fCalibFDC1Track;
  if (fCalibFDC2Hit)    delete fCalibFDC2Hit;
  if (fCalibFDC2Track)    delete fCalibFDC2Track;
}

void AnaFDC::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("fdc1X",&fdc1X,"fdc1X/D");
  tree->Branch("fdc1A",&fdc1A,"fdc1A/D");
  tree->Branch("fdc1Y",&fdc1Y,"fdc1Y/D");
  tree->Branch("fdc1B",&fdc1B,"fdc1B/D");
  tree->Branch("fdc1C",&fdc1C,"fdc1C/D");
  tree->Branch("fdc2X",&fdc2X,"fdc2X/D");
  tree->Branch("fdc2A",&fdc2A,"fdc2A/D");
  tree->Branch("fdc2Y",&fdc2Y,"fdc2Y/D");
  tree->Branch("fdc2B",&fdc2B,"fdc2B/D");
  tree->Branch("fdc2C",&fdc2C,"fdc2C/D");

  tree->Branch("fdc1ex",&fdc1ex,"fdc1ex/D");
  tree->Branch("fdc1ey",&fdc1ey,"fdc1ey/D");
  tree->Branch("fdc1ez",&fdc1ez,"fdc1ez/D");
}


void AnaFDC::LoadTDCDistribution(const char* filename)
{
  char myname[128];
  TFile *fdcin = new TFile(filename,"READ"); 
  gROOT->cd();
  TH1F *hist = NULL;

  for(int i=0;i<7;i++){
    sprintf(myname,"fdc1_ftdc_corr_%d",i);

    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionX();
    fCalibFDC1Track->SetTDCDistribution(hist,i*2);
    delete hist; hist = NULL;
    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionY();
    fCalibFDC1Track->SetTDCDistribution(hist,i*2+1);
    delete hist; hist = NULL;

    sprintf(myname,"fdc2_ftdc_corr_%d",i);

    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionX();
    fCalibFDC2Track->SetTDCDistribution(hist,i*2);
    delete hist; hist = NULL;
    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionY();
    fCalibFDC2Track->SetTDCDistribution(hist,i*2+1);
    delete hist; hist = NULL;
  }
}
