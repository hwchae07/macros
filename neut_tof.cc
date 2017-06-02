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
  
  TH2 *hist1 = new TH2D("hist1","hist1",400,0.5,400.5,1000,-890,-840);
  TH2 *hsum1 = new TH2D("hsum1","hsum1",400,0.5,400.5,1000,-890,-840);


  for(Int_t run=146;run<=154;run++)
    {
      //run = 146;
      TFile *file = new TFile(Form("./root/run%04d.root.BeamPla",run),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("RunNum","EventNum");
      tree->AddFriend("NeuLAND",Form("./root/run%04d.root.NeuLAND",run));
      tree->GetFriend("NeuLAND")->BuildIndex("RunNum","EventNum");

      //beam cut 22Ne, 21F and fragment cut 22Ne//
      TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
      TCutG *b22ne, *b21f, *f22ne;
      f_cut->GetObject("b22ne",b22ne);
      f_cut->GetObject("b21f",b21f);
      f_cut->GetObject("f22ne",f22ne);
      //beam cut 22Ne, 21F and fragment cut 22Ne//


      tree->SetAlias("FL_tgt_neuland","sqrt(neulandX**2 + neulandY**2 + neulandZ**2)");
      
      tree->Draw("neulandTA-T13-FL_tgt_neuland/299.792:neulandID>>hist1(400,0.5,400.5,1000,-890,-840)","","goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);

    }
 
  return hsum1;
}

void savehist()
{
  TH2 *h1 = make2dhist();
  TFile *file = new TFile("./root/neulandTOF.root","recreate");
  h1->Write();
  file->Close();
}


void gammafit()
{
  TFile *file = new TFile("./root/neulandTOF.root","r");
  TH2 *h1 = (TH2D*)file->Get("hsum1");
  TH1 *hist;
  TSpectrum *ts = new TSpectrum();
  TF1 *fit1;
  TF1 *fit2 = new TF1("gausP2",gausP2,-313,-311,5);

  TCanvas *c1 = new TCanvas("c1","c1",500,500);

  h1->Draw("colz");



}


