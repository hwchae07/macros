#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TSpectrum.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include <iostream>
#include <fstream>

void tcal_CATANA()
{
  TFile *file = new TFile("./root/run0319.root.CATANA");
  TTree *tree = (TTree*)file->Get("CATANA");

  TH1 *ht = new TH1D("ht","ht",8000,0,8000);
  TSpectrum *ts1 = new TSpectrum();

  TString detName = "CATANA";
  
  ofstream file1(Form("./parameter/tcal_%s.dat",detName.Data()));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  c1->cd(1);
  tree->Draw("catanaRawTDC-catanaRawTRef>>ht");
  ts1->Search(ht,2,"",0.005);
  cout<<ts1->GetNPeaks()<<endl;
  Double_t *peakX1 = new Double_t[ts1->GetNPeaks()];
  Int_t *index1 = new Int_t[ts1->GetNPeaks()];
  Double_t *ch1 = new Double_t[ts1->GetNPeaks()];
  Double_t *err1 = new Double_t[ts1->GetNPeaks()];
  Double_t *time1 = new Double_t[ts1->GetNPeaks()];//psec
  TMath::Sort(ts1->GetNPeaks(),ts1->GetPositionX(),index1,false);
  for(Int_t i=0;i<ts1->GetNPeaks();i++)
    {
      peakX1[i] = (ts1->GetPositionX())[index1[i]];
      time1[i] = 10 * (i+1);
      //cout<<i<<"\t"<<index1[i]<<"\t"<<peakX1[i]<<endl;
    }

  c1->Print(Form("./fig/RawTDC_%s.pdf",detName.Data()));
  //c1->Close();

  
   
  TH1 *hist;
  TF1 *fit1, *fit2;
  Double_t min,max,center;
  Double_t range = 30;
  Double_t width = 5;
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,900);
  c2->Divide(3,3);
  c2->Print(Form("./fig/traw_%s.pdf[",detName.Data()));
  for(Int_t id=0;id<ts1->GetNPeaks();id++)
    {
      c2->cd((id%9)+1);
      center = peakX1[id];
      max = center + range;
      min = center - range;
      bin = range;
      
      tree->Draw(Form("catanaRawTDC-catanaRawTRef>>ht%d(%lf,%lf,%lf)",id+1,bin,min,max));
      hist = (TH1D*)gDirectory->Get(Form("ht%d",id+1));
      hist->SetTitle(Form("Peak %d",id+1));
      hist->GetXaxis()->SetTitle("channel (ch)");
      hist->Fit("gaus","rq","",center-width,center+width);
      fit1 = (TF1*)hist->GetFunction("gaus");
      ch1[id] = fit1->GetParameter(1);
      err1[id] = fit1->GetParameter(2);
      cout<<"Peak "<<id+1<<" is fitting..."<<endl;
      c2->Update();

      if((id+1)%9==0)
	c2->Print(Form("./fig/traw_%s.pdf",detName.Data()));
    }
  c2->Print(Form("./fig/traw_%s.pdf]",detName.Data()));


  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  TGraphErrors *gr1 = new TGraphErrors(ts1->GetNPeaks(),ch1,time1,err1,0);
  gr1->Fit("pol1","q");
  gr1->SetMarkerStyle(24);
  gr1->SetTitle(Form("%s Raw TDC ",detName.Data()));
  gr1->GetXaxis()->SetTitle("channel (ch)");
  gr1->GetYaxis()->SetTitle("time (nsec)");
  gr1->Draw("AP");
  c4->Update();
  c4->Print(Form("./fig/tcal_%s.pdf",detName.Data()));

  TF1 *func1 = (TF1*)gr1->GetFunction("pol1");
  file1<<"Left"<<endl;
  file1<<"offset\t\t"<<func1->GetParameter(0)<<std::endl;
  file1<<"slope\t\t"<<func1->GetParameter(1)<<std::endl;
  file1<<"chi2/NDF\t"<<func1->GetChisquare()/func1->GetNDF()<<std::endl;
 
  
  
  delete [] peakX1,index1,ch1,err1;
  file1.close();


}
