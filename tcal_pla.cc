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

void tcal_pla()
{
  TFile *file = new TFile("./root/run0302.root.BeamPla");
  TTree *tree = (TTree*)file->Get("BeamPla");

  TH1 *htl = new TH1D("htl","htl",33000,0,66000);
  TH1 *htr = new TH1D("htr","htr",33000,0,66000);


  Int_t detNum;
  TString detNamelist[5] = {"f3Pla","f5Pla","f7Pla","f13Pla1","f13Pla2"};
  TString detName;
  
  cout<<"Choose detector"<<endl;
  cout<<"F3  Plastic  : 0"<<endl;
  cout<<"F5  Plastic  : 1"<<endl;
  cout<<"F7  Plastic  : 2"<<endl;
  cout<<"F13 Plastic1 : 3"<<endl;
  cout<<"F13 Plastic2 : 4"<<endl;
 
  cin>>detNum;

  detName = detNamelist[detNum];

  ofstream file1(Form("./parameter/tcal_%s.dat",detName.Data()));
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);

  c1->cd(1);
  tree->Draw(Form("%sTLRaw>>htl",detName.Data()));
  c1->cd(2);
  tree->Draw(Form("%sTRRaw>>htr",detName.Data()));

  Int_t peakN1 = 0;
  Int_t peakN2 = 0;
  Double_t peakX1[200]={0,};
  Double_t peakX2[200]={0,};


  for(Int_t i=0;i<=htl->GetNbinsX();i++)
    {
      if(htl->GetBinContent(i)<=30 && htl->GetBinContent(i+1)>30)
	{
	  peakX1[peakN1] = htl->GetBinCenter(i);
	  peakN1++;
	  //std::cout<<peakN<< " : " <<htl->GetBinCenter(i)<<endl;
	  i+=100;
	}
    }

  
  for(Int_t i=0;i<=htr->GetNbinsX();i++)
    {
      if(htr->GetBinContent(i)<=30 && htr->GetBinContent(i+1)>30)
	{
	  peakX2[peakN2] = htr->GetBinCenter(i);
	  peakN2++;
	  //std::cout<<peakN<< " : " <<htl->GetBinCenter(i)<<endl;
	  i+=100;
	}
    }

  std::cout<<"Left Raw TDC Peak number is "<<peakN1<<endl;
  std::cout<<"Right Raw TDC Peak number is "<<peakN2<<endl;
  
  Double_t *ch1 = new Double_t[peakN1];
  Double_t *err1 = new Double_t[peakN1];
  Double_t *time1 = new Double_t[peakN1];
  Double_t *res1 = new Double_t[peakN1];

  Double_t *ch2 = new Double_t[peakN2];
  Double_t *err2 = new Double_t[peakN2];
  Double_t *time2 = new Double_t[peakN2];
  Double_t *res2 = new Double_t[peakN2];


  TCanvas *c2 = new TCanvas("c2","c2",1600,960);
  c2->Divide(5,4);
  c2->Print(Form("./fig/tlraw_%s.pdf[",detName.Data()));
  
  TH1 *hist;
  TF1 *fit1;
  Double_t min,max,center,bin;
  Double_t range = 30;
  Double_t width = 10;

  
  for(Int_t id=0;id<peakN1;id++)
    {
      time1[id] = (id+1) * 10;
      c2->cd((id%20)+1);
      center = peakX1[id];
      max = center + range;
      min = center - range;
      bin = range;
      //cout<<center[id]<<endl;
      tree->Draw(Form("%sTLRaw>>htl%d(%lf,%lf,%lf)",detName.Data(),id+1,bin,min,max),Form("abs(%sTLRaw-%lf)<%lf",detName.Data(),center,range));
      hist = (TH1D*)gDirectory->Get(Form("htl%d",id+1));
      hist->SetTitle(Form("Peak %d",id+1));
      hist->GetXaxis()->SetTitle("channel (ch)");
      Double_t ndf = hist->GetEntries();

      
      hist->Fit("gaus","rq","",center-width,center+width);
      fit1 = (TF1*)hist->GetFunction("gaus");
      ch1[id] = fit1->GetParameter(1);
      err1[id] = fit1->GetParameter(2)/sqrt(ndf);

      for(Int_t j=0;j<3;j++)
	{
	  fit1->SetParameter(1,ch1[id]);
	  fit1->SetParameter(2,err1[id]);
	  hist->Fit(fit1,"rq","",ch1[id]-width,ch1[id]+width);
	  ch1[id] = fit1->GetParameter(1);
	  err1[id] = fit1->GetParameter(2)/sqrt(ndf);
	}
      cout<<"Peak "<<id+1<<" is fitting..."<<endl;
      c2->Update();

      if((id+1)%20==0||id==peakN1-1)
	{
	  c2->Print(Form("./fig/tlraw_%s.pdf",detName.Data()));
	  c2->Clear();
	  c2->Divide(5,4);
	}
      
    }
  c2->Print(Form("./fig/tlraw_%s.pdf]",detName.Data()));



  
  TCanvas *c3 = new TCanvas("c3","c3",1600,960);
  c3->Divide(5,4);
  c3->Print(Form("./fig/trraw_%s.pdf[",detName.Data()));
  
  for(Int_t id=0;id<peakN2;id++)
    {
      time2[id] = (id+1) * 10;
      c3->cd((id%20)+1);
      center = peakX2[id];
      max = center + range;
      min = center - range;
      bin = range;
      //cout<<center[id]<<endl;
      tree->Draw(Form("%sTRRaw>>htr%d(%lf,%lf,%lf)",detName.Data(),id+1,bin,min,max),Form("abs(%sTRRaw-%lf)<%lf",detName.Data(),center,range));
      hist = (TH1D*)gDirectory->Get(Form("htr%d",id+1));
      hist->SetTitle(Form("Peak %d",id+1));
      hist->GetXaxis()->SetTitle("channel (ch)");
      Double_t ndf = hist->GetEntries();
      
      hist->Fit("gaus","rq","",center-width,center+width);
      fit1 = (TF1*)hist->GetFunction("gaus");
      ch2[id] = fit1->GetParameter(1);
      err2[id] = fit1->GetParameter(2)/sqrt(ndf);

      for(Int_t j=0;j<3;j++)
	{
	  fit1->SetParameter(1,ch2[id]);
	  fit1->SetParameter(2,err2[id]);
	  hist->Fit(fit1,"rq","",ch2[id]-width,ch2[id]+width);
	  ch2[id] = fit1->GetParameter(1);
	  err2[id] = fit1->GetParameter(2)/sqrt(ndf);
	}
      cout<<"Peak "<<id+1<<" is fitting..."<<endl;
      c3->Update();

      if((id+1)%20==0||id==peakN2-1)
	{
	  c3->Print(Form("./fig/trraw_%s.pdf",detName.Data()));
	  c3->Clear();
	  c3->Divide(5,4);
	}
      
    }
  c3->Print(Form("./fig/trraw_%s.pdf]",detName.Data()));
  



  

  TCanvas *c4 = new TCanvas("c4","c4",800,600);
  c4->Print(Form("./fig/tcal_%s.pdf[",detName.Data()));
  TGraphErrors *gr1 = new TGraphErrors(peakN1,ch1,time1,err1,0);
  gr1->Fit("pol1","q");
  gr1->SetMarkerStyle(24);
  gr1->SetTitle(Form("%s Raw TimeL ",detName.Data()));
  gr1->GetXaxis()->SetTitle("channel (ch)");
  gr1->GetYaxis()->SetTitle("time (nsec)");
  gr1->Draw("AP");
  c4->Update();
  c4->Print(Form("./fig/tcal_%s.pdf",detName.Data()));

  TGraphErrors *gr2 = new TGraphErrors(peakN2,ch2,time2,err2,0);
  gr2->Fit("pol1","q");
  gr2->SetMarkerStyle(24);
  gr2->SetTitle(Form("%s Raw TimeR ",detName.Data()));
  gr2->GetXaxis()->SetTitle("channel (ch)");
  gr2->GetYaxis()->SetTitle("time (nsec)");
  gr2->Draw("AP");
  c4->Update();
  c4->Print(Form("./fig/tcal_%s.pdf",detName.Data()));

  c4->Print(Form("./fig/tcal_%s.pdf]",detName.Data()));


  


  //file1<<detName.Data()<<endl;

  TCanvas *c5 = new TCanvas("c5","c5",800,600);

  c5->Print(Form("./fig/tdc_res_%s.pdf[",detName.Data()));
  TF1 *func1 = (TF1*)gr1->GetFunction("pol1");
  for(Int_t i=0;i<peakN1;i++)
    res1[i] = time1[i] - func1->Eval(ch1[i]);

  TGraph *gr3 = new TGraph(peakN1,ch1,res1);
  gr3->Draw("AL");
  gr3->SetTitle(Form("%s Residual Time (L);channel (ch);Time (ns)",detName.Data()));
  c5->Update();
  c5->Print(Form("./fig/tdc_res_%s.pdf",detName.Data()));
  
  file1<<"Left"<<endl;
  file1<<"offset\t\t"<<func1->GetParameter(0)<<std::endl;
  file1<<"slope\t\t"<<func1->GetParameter(1)<<std::endl;
  file1<<"chi2/NDF\t"<<func1->GetChisquare()/func1->GetNDF()<<std::endl;

  TF1 *func2 = (TF1*)gr2->GetFunction("pol1");
  for(Int_t i=0;i<peakN2;i++)
    res2[i] = time2[i] - func2->Eval(ch2[i]);

  TGraph *gr4 = new TGraph(peakN2,ch2,res2);
  gr4->Draw("AL");
  gr4->SetTitle(Form("%s Residual Time (R);channel (ch);Time (ns)",detName.Data()));
  c5->Update();
  c5->Print(Form("./fig/tdc_res_%s.pdf",detName.Data()));

  c5->Print(Form("./fig/tdc_res_%s.pdf]",detName.Data()));


  
  //histogram for residue...
  TCanvas *c6 = new TCanvas("c6","c6",800,600);
  c6->Print(Form("./fig/tdc_res_hist_%s.pdf[",detName.Data()));
  TH1 *h1 = new TH1D("h1","h1",15,-0.03,0.03);
  TH1 *h2 = new TH1D("h2","h2",15,-0.03,0.03);
  for(Int_t i=0;i<peakN1;i++)
    h1->Fill(res1[i]);
  for(Int_t i=0;i<peakN2;i++)
    h2->Fill(res2[i]);
  h1->SetTitle(Form("%s Residual Time (L);Residual Time (ns)",detName.Data()));
  h2->SetTitle(Form("%s Residual Time (R);Residual Time (ns)",detName.Data()));
  h1->Draw();
  c6->Update();
  c6->Print(Form("./fig/tdc_res_hist_%s.pdf",detName.Data()));
  h2->Draw();
  c6->Update();
  c6->Print(Form("./fig/tdc_res_hist_%s.pdf",detName.Data()));

  c6->Print(Form("./fig/tdc_res_hist_%s.pdf]",detName.Data()));

  
  
  file1<<"Right"<<endl;
  file1<<"offset\t\t"<<func2->GetParameter(0)<<std::endl;
  file1<<"slope\t\t"<<func2->GetParameter(1)<<std::endl;
  file1<<"chi2/NDF\t"<<func2->GetChisquare()/func2->GetNDF()<<std::endl;

  
  
  delete [] ch1,err1,time1,ch2,err2,time2;
  file1.close();

  c1->Close();
  c2->Close();
  c3->Close();
  c4->Close();
  c4->Close();
  c5->Close();
  c6->Close();

}
