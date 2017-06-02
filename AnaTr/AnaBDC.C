#include "AnaBDC.H"

#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "TArtDCTrack.hh"

AnaBDC::AnaBDC():
  AnaModule("BDC"),
  fBDCParameters(NULL),
  fCalibBDC1Hit(NULL),
  fCalibBDC1Track(NULL),
  fCalibBDC2Hit(NULL),
  fCalibBDC2Track(NULL) {
  ;}

AnaBDC::~AnaBDC(){    
  DeleteAll();}

void AnaBDC::InitParameter(){
  fBDCParameters = TArtSAMURAIParameters::Instance();
  fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC1.xml");
  fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC2.xml");
  parLoaded = true;}

void AnaBDC::InitDetector(){
  if (!parLoaded) return;
  fCalibBDC1Hit = new TArtCalibBDC1Hit;
  fCalibBDC1Track = new TArtCalibBDC1Track;
  fCalibBDC2Hit = new TArtCalibBDC2Hit;
  fCalibBDC2Track = new TArtCalibBDC2Track;

  LoadTDCDistribution((char*)"db/dc/s027_run0167.root");
  detLoaded = true;}

void AnaBDC::Analysis(){
  anaFlag = true;
  if (!detLoaded) return;

  fCalibBDC1Hit->ClearData();
  fCalibBDC1Track->ClearData();
  fCalibBDC2Hit->ClearData();
  fCalibBDC2Track->ClearData();

  fCalibBDC1Hit->ReconstructData();
  fCalibBDC1Track->ReconstructData();
  fCalibBDC2Hit->ReconstructData();
  fCalibBDC2Track->ReconstructData();

  if (fCalibBDC1Track->GetNumDCTrack() < 2) {
    anaFlag = false; return; }
  if (fCalibBDC2Track->GetNumDCTrack() < 2) {
    anaFlag = false; return; }

  Double_t chi2X = 100000;
  Double_t chi2Y = 100000;
  bdc1X = -9999; bdc1Y = -9999;

  for (Int_t i = 0 ; i < fCalibBDC1Track->GetNumDCTrack() && (bdc1X == -9999 || bdc1Y == -9999.) ; i++){
    TArtDCTrack *trk = fCalibBDC1Track->GetDCTrack(i);
    Double_t chi2 = trk->GetChi2();
    Int_t ndf = trk->GetNDF();
    Double_t posx = trk->GetPosition(0);
    Double_t posy = trk->GetPosition(1);
    Double_t angx = trk->GetAngle(0);
    Double_t angy = trk->GetAngle(1);
    if (posx < -5000){
      if (chi2/(Double_t)ndf < chi2Y){
	chi2Y = chi2/(Double_t)ndf;
	bdc1Y = posy; bdc1B = angy; bdc1CY = chi2Y;}}
    else if (posy < -5000){
      if (chi2/(Double_t)ndf < chi2X){
	chi2X = chi2/(Double_t)ndf;
	bdc1X = posx; bdc1A = angx; bdc1CX = chi2X;}}}
  //  bdc1C = sqrt(chi2X*chi2X + chi2Y*chi2Y);

  if (bdc1X == -9999 || bdc1Y == -9999) {
    anaFlag = false; return; }

  chi2X = 100000;
  chi2Y = 100000;
  bdc2X = -9999; bdc2Y = -9999;

  for (Int_t i = 0 ; i < fCalibBDC2Track->GetNumDCTrack() && (bdc2X == -9999 || bdc2Y == -9999.) ; i++){
    TArtDCTrack *trk = fCalibBDC2Track->GetDCTrack(i);
    Double_t chi2 = trk->GetChi2();
    Int_t ndf = trk->GetNDF();
    Double_t posx = trk->GetPosition(0);
    Double_t posy = trk->GetPosition(1);
    Double_t angx = trk->GetAngle(0);
    Double_t angy = trk->GetAngle(1);
    if (posx < -5000){
      if (chi2/(Double_t)ndf < chi2Y){
	chi2Y = chi2/(Double_t)ndf;
	bdc2Y = posy; bdc2B = angy; bdc2CY = chi2Y;}}
    else if (posy < -5000){
      if (chi2/(Double_t)ndf < chi2X){
	chi2X = chi2/(Double_t)ndf;
	bdc2X = posx; bdc2A = angx; bdc2CX = chi2X;}}}
  //  bdc2C = sqrt(chi2X*chi2X + chi2Y*chi2Y);

  if (bdc2X == -9999 || bdc2Y == -9999) {
    anaFlag = false; return; }
}

void AnaBDC::DeleteAll(){
  if (fCalibBDC1Hit)    delete fCalibBDC1Hit;
  if (fCalibBDC1Track)    delete fCalibBDC1Track;
  if (fCalibBDC2Hit)    delete fCalibBDC2Hit;
  if (fCalibBDC2Track)    delete fCalibBDC2Track;
}

void AnaBDC::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("bdc1X",&bdc1X,"bdc1X/D");
  tree->Branch("bdc1A",&bdc1A,"bdc1A/D");
  tree->Branch("bdc1Y",&bdc1Y,"bdc1Y/D");
  tree->Branch("bdc1B",&bdc1B,"bdc1B/D");
  tree->Branch("bdc1CX",&bdc1CX,"bdc1CX/D");
  tree->Branch("bdc1CY",&bdc1CY,"bdc1CY/D");
  tree->Branch("bdc2X",&bdc2X,"bdc2X/D");
  tree->Branch("bdc2A",&bdc2A,"bdc2A/D");
  tree->Branch("bdc2Y",&bdc2Y,"bdc2Y/D");
  tree->Branch("bdc2B",&bdc2B,"bdc2B/D");
  tree->Branch("bdc2CX",&bdc2CX,"bdc2CX/D");
  tree->Branch("bdc2CY",&bdc2CY,"bdc2CY/D");

  tree->SetAlias("tgtX","bdc1X + 1853. / 1000. * (bdc2X - bdc1X)");
  tree->SetAlias("tgtY","bdc1Y + 1853. / 1000. * (bdc2Y - bdc1Y)");
  tree->SetAlias("tgtA","tan( (bdc2X - bdc1X) / 1000. )");
  tree->SetAlias("tgtB","tan( (bdc2Y - bdc1Y) / 1000. )");
}


void AnaBDC::LoadTDCDistribution(const char* filename)
{
  char myname[128];
  TFile *fdcin = new TFile(filename,"READ"); 
  gROOT->cd();
  TH1F *hist = NULL;

  for (Int_t i = 0 ; i < 4 ; i++){
    sprintf(myname,"bdc1_ftdc_corr_%d",i);
    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionX();
    fCalibBDC1Track->SetTDCDistribution(hist,i*2);
    delete hist; hist = NULL;
    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionY();
    fCalibBDC1Track->SetTDCDistribution(hist,i*2+1);
    delete hist; hist = NULL;

    sprintf(myname,"bdc2_ftdc_corr_%d",i);
    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionX();
    fCalibBDC2Track->SetTDCDistribution(hist,i*2);
    delete hist; hist = NULL;
    hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionY();
    fCalibBDC2Track->SetTDCDistribution(hist,i*2+1);
    delete hist; hist = NULL;
  }
}
