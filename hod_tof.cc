#include "TFile.h"
#include "TTree.h"


//void hod_tof()
{

  TH2 *hist1 = new TH2D("hist1","hist1",200,80,90,200,200,550);
  TH2 *hsum1 = new TH2D("hsum1","HOD Time VS Charge for {}^{32}Ne;Time;Charge",200,80,90,200,200,550);
  TH2 *hist2 = new TH2D("hist2","hist2",200,-60,60,200,-60,60);
  TH2 *hsum2 = new TH2D("hsum2","FDC1 X VS Y for {}^{32}Ne;X (mm);Y (mm)",200,-60,60,200,-60,60);
  TH2 *hist3 = new TH2D("hist3","hist3",200,0,300,200,-300,300);
  TH2 *hsum3 = new TH2D("hsum3","FDC2 X VS Y for {}^{32}Ne;X (mm);Y (mm)",200,0,300,200,-300,300);  
  TH2 *hist4 = new TH2D("hist4","hist4",24,0.5,24.5,50,83,89);
  TH2 *hsum4 = new TH2D("hsum4","HOD ID VS #DeltaT for {}^{32}Ne;ID;#DeltaT",24,0.5,24.5,50,83,89);

  TH1 *hist5 = new TH1D("hist5","hist5",200,-228,-222);
  TH1 *hsum5 = new TH1D("hsum5","#DeltaT for {}^{32}Ne;#DeltaT",200,-228,-222);
  TCanvas *c1 = new TCanvas("c1","c1",1600,400);
  c1->Divide(4,1);
  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  
  Int_t runStart = 275;
  Int_t runEnd = 275;
  //for(Int_t runNum = runStart ; runNum<=runEnd ; runNum++)
  //{
      runNum = 275;
      TString filename = Form("run%04d.root",runNum);
  
      TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("EventNum");
      tree->AddFriend("BDC",Form("./root/%s.BDC",filename.Data()));
      tree->GetFriend("BDC")->BuildIndex("EventNum");
      tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
      tree->GetFriend("FDC")->BuildIndex("EventNum");
      tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
      tree->GetFriend("HOD")->BuildIndex("EventNum");
      tree->AddFriend("Brho",Form("./root/%s.Brho",filename.Data()));
      tree->GetFriend("Brho")->BuildIndex("EventNum");
  
      //Beam PID//
      tree->SetAlias("tofc","TOF3_13-3*dT5");
      tree->SetAlias("elossc","icbEloss-0.1*tofc-40");
      TCut cut_na = "TMath::Abs(elossc-22.5)<2.5";
      TCut cut_34na = cut_na&&"TMath::Abs(tofc+388)<3";
      TCut cut_ne = "abs(elossc-17.8)<2.2";
      TCut cut_32ne = cut_ne&&"abs(tofc+379)<3";
      //Beam PID//

      //Target size//
      Double_t bdcWidth = 120.;
      Double_t bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
      Double_t bdc1Tgt = bdc2Tgt + 1000;
      tree->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
      tree->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
      TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";  
      //Target size//

      //Fragment PID//
      TCut cut_f32ne = "abs(hodTA-85.6)<0.8&&abs(hodQA-410)<10";
      //Fragment PID//

      //FDC1 & FDC2//
      TCut cut_fdc1 = "sqrt(fdc1X**2 + fdc1Y**2) < 40";
      TCut cut_fdc2 = "((fdc2X-150)/100)**2 + (fdc2Y/200)**2 < 1";
      //FDC1 & FDC2//
      TCut cut_total = cut_32ne&&cut_reaction&&cut_f32ne;

      c1->cd(1);
      tree->Draw("hodQA:hodTA>>hist1(200,80,90,200,200,550)",cut_total);
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      hsum1->Draw("colz");

      c1->cd(2);
      tree->Draw("fdc1Y:fdc1X>>hist2(200,-60,60,200,-60,60)",cut_total&&cut_fdc1);
      gDirectory->GetObject("hist2",hist2);
      hsum2->Add(hist2);
      hsum2->Draw("colz");

      c1->cd(3);
      tree->Draw("fdc2Y:fdc2X>>hist3(200,0,300,200,-300,300)",cut_total&&cut_fdc1&&cut_fdc2);
      gDirectory->GetObject("hist3",hist3);
      hsum3->Add(hist3);
      hsum3->Draw("colz");
      
      c1->cd(4);
      tree->Draw("hodTA:hodID>>hist4(24,0.5,24.5,50,83,89)",cut_total&&cut_fdc1&&cut_fdc2);
      gDirectory->GetObject("hist4",hist4);
      hsum4->Add(hist4);
      hsum4->Draw("colz");
      //}

      c2->cd();
      tree->Draw("hodTA-T13>>hist5(200,-228,-222)",cut_total&&cut_fdc1&&cut_fdc2&&"hodID==11");
      gDirectory->GetObject("hist5",hist5);
      hsum5->Add(hist5);
      hsum5->Draw();
      hsum5->Fit("gaus","r","",-225,-224);
}
