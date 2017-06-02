#include "TFile.h"
#include "TTree.h"
#include "TCutG.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include <iostream>
#include <fstream>

Double_t gausP2(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t result;

  result = par[0] * exp(-0.5*((xx-par[1])/par[2])*(xx-par[1])/par[2]);
  if(xx>par[1])
    result += par[3] * (xx-par[1]) * exp(-par[4]*(xx-par[1]));
  result += par[5];
  return result;
}

TH2* make2dhist(){
  Double_t fl_nebula;
  
  TH2 *hist1 = new TH2D("hist1","hist1",144,0.5,144.5,50,-313.5,-306);
  TH2 *hsum1 = new TH2D("hsum1","hsum1",144,0.5,144.5,50,-313.5,-306);

  for(Int_t run=146;run<=154;run++)
    {
      //run = 146;
      TFile *file = new TFile(Form("./root/run%04d.root.BeamPla",run),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("RunNum","EventNum");
      tree->AddFriend("NEBULA",Form("./root/run%04d.root.NEBULA",run));
      tree->GetFriend("NEBULA")->BuildIndex("RunNum","EventNum");

      //beam cut 22Ne, 21F and fragment cut 22Ne//
      TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
      TCutG *b22ne, *b21f, *f22ne;
      f_cut->GetObject("b22ne",b22ne);
      f_cut->GetObject("b21f",b21f);
      f_cut->GetObject("f22ne",f22ne);
      //beam cut 22Ne, 21F and fragment cut 22Ne//


      tree->SetAlias("FL_tgt_nebula","sqrt(nebulaX**2 + nebulaY**2 + nebulaZ**2)");
      
      tree->Draw("nebulaTA-T13-FL_tgt_nebula/299.792:nebulaID>>hist1(144,0.5,144.5,50,-313.5,-306)","","goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);

    }
 
  return hsum1;
}

void savehist()
{
  TH2 *h1 = make2dhist();
  TFile *file = new TFile("./root/nebulaTOF.root","recreate");
  h1->Write();
  file->Close();
}


void gammafit()
{
  TFile *file = new TFile("./root/nebulaTOF.root","r");
  TH2 *h1 = (TH2D*)file->Get("hsum1");
  TH1 *hist;
  TSpectrum *ts = new TSpectrum();
  TF1 *fit1;
  TF1 *fit2 = new TF1("gausP2",gausP2,-313,-311,5);

  TCanvas *c1 = new TCanvas("c1","c1",1500,900);
  c1->Divide(5,3);

  Double_t tof_c, tof_h, tof_l;
  Double_t tof_r = 1;
  Double_t range_min,range_max;
  Int_t canvasNum = 1;

  ofstream fout("./dat/nebula_gamma_offset.dat");
  
  for(Int_t id=1;id<=120;id++)
    {
      
      c1->cd((id-1)%15+1);
      h1->ProjectionY(Form("h%02d",id),id,id);
      gDirectory->GetObject(Form("h%02d",id),hist);
      hist->SetTitle(Form("ID %02d;TOF",id));
      ts->Search(hist,1,"",0.5);
      tof_c = TMath::MinElement(ts->GetNPeaks(),ts->GetPositionX());
      tof_l = tof_c - tof_r;
      tof_h = tof_c + tof_r;

      hist->Fit("gaus","RQ","",tof_l,tof_h);
      fit1 = hist->GetFunction("gaus");
      /*
      fit2->SetParameter(0,fit1->GetParameter(0));
      fit2->SetParameter(1,fit1->GetParameter(1));
      fit2->SetParameter(2,fit1->GetParameter(2));
      fit2->SetParameter(3,10);
      fit2->SetParameter(4,1);

      hist->Fit("gausP2","RQ","",tof_l,tof_h);
      */

      fout<<fit1->GetParameter(1)<<std::endl;
      
      hist->Draw();
      c1->Update();
      if(!(id%15))
	getchar();

    }

  fout.close();
}


void check2dhist(){
  Double_t fl_nebula;
  
  TH2 *hist1 = new TH2D("hist1","hist1",120,0.5,120.5,500,-20,80);
  TH2 *hsum1 = new TH2D("hsum1","hsum1",120,0.5,120.5,500,-20,80);

  TH2 *hist2 = new TH2D("hist2","hist2",120,0.5,120.5,500,-280,-180);
  TH2 *hsum2 = new TH2D("hsum2","hsum2",120,0.5,120.5,500,-280,-180);

  ifstream fin("./dat/nebula_gamma_offset.dat");
  Double_t nebula_offset[120];
  for(Int_t i=0;i<120;i++)
    fin>>nebula_offset[i];
  
  TString nebula_tof = "0";

  for(Int_t id=1;id<=120;id++)
    {
      nebula_tof += "+(nebulaID==";
      nebula_tof += id;
      nebula_tof += ")*(";
      nebula_tof += -nebula_offset[id-1];
      nebula_tof += ")";
    }
  //  std::cout<<nebula_tof<<std::endl;

  TCanvas *c1 = new TCanvas("c1","c1",500,500);
  
  for(Int_t run=146;run<=154;run++)
    {
      TFile *file = new TFile(Form("./root/run%04d.root.BeamPla",run),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("RunNum","EventNum");
      tree->AddFriend("NEBULA",Form("./root/run%04d.root.NEBULA",run));
      tree->GetFriend("NEBULA")->BuildIndex("RunNum","EventNum");

      //beam cut 22Ne, 21F and fragment cut 22Ne//
      TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
      TCutG *b22ne, *b21f, *f22ne;
      f_cut->GetObject("b22ne",b22ne);
      f_cut->GetObject("b21f",b21f);
      f_cut->GetObject("f22ne",f22ne);
      //beam cut 22Ne, 21F and fragment cut 22Ne//


      tree->SetAlias("FL_tgt_nebula","sqrt(nebulaX**2 + nebulaY**2 + nebulaZ**2)");

      /*
	tree->Draw("nebulaTA-T13:nebulaID>>hist2(120,0.5,120.5,500,-280,-180)","","colz");
	gDirectory->GetObject("hist2",hist2);
	hsum2->Add(hist2);
	hsum2->Draw("colz");
      */
 
      tree->Draw(Form("nebulaTA-T13-FL_tgt_nebula/299.792+%s:nebulaID>>hist1(120,0.5,120.5,500,-20,80)",nebula_tof.Data()),"","colz");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      hsum1->Draw("colz");
 
      c1->Update();

    }
    

  
}
