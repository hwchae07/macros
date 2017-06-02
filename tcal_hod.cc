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

void tcal_hod()
{
  TFile *file = new TFile("./root/run0302.root.HOD");
  TTree *tree = (TTree*)file->Get("HOD");

  TH1 *ht1 = new TH1D("ht1","ht1",4950,51.5,5001.5);
  TH1 *ht2 = new TH1D("ht2","ht2",4950,51.5,5001.5);
  TSpectrum *ts1 = new TSpectrum();
  TSpectrum *ts2 = new TSpectrum();

  TString detName = "HOD";

  ofstream fout("./parameter/tcal_HOD.dat");
  TF1 *fit;
  fout<<setw(10)<<"up_off\t"<<setw(10)<<"up_slope\t"<<setw(10)<<"down_off\t"<<setw(10)<<"down_slope"<<std::endl;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  TCanvas *c3 = new TCanvas("c3","c3",1200,600);
  c3->Divide(2,1);
  
  c1->Print("./fig/RawTDC_HOD.pdf[");
  c2->Print("./fig/tcal_HOD.pdf[");
  c3->Print("./fig/tcal_HOD_res.pdf[");
  //Int_t id=1;

  for(Int_t id=1;id<25;id++)
    {
  
      c1->cd(1);
      tree->Draw("hodTURaw>>ht1(4950,51.5,5001.5)",Form("hodID==%d",id));
      ht1 = (TH1D*)gDirectory->Get("ht1");
      ht1->SetNameTitle(Form("HOD Raw_up TDC ID%d",id),Form("HOD Raw_up TDC ID%d",id));
      ts1->Search(ht1,2,"",0.1);
      cout<<"ID"<<id<<"(u): "<<ts1->GetNPeaks()<<endl;
      Double_t *peakX1 = new Double_t[ts1->GetNPeaks()];
      Int_t *index1 = new Int_t[ts1->GetNPeaks()];
      Double_t *ch1 = new Double_t[ts1->GetNPeaks()];
      Double_t *err1 = new Double_t[ts1->GetNPeaks()];
      Double_t *time1 = new Double_t[ts1->GetNPeaks()];//psec
      Double_t *res1 = new Double_t[ts1->GetNPeaks()];
      TMath::Sort(ts1->GetNPeaks(),ts1->GetPositionX(),index1,false);
      for(Int_t i=0;i<ts1->GetNPeaks();i++)
	{
	  peakX1[i] = (ts1->GetPositionX())[index1[i]];
	  time1[i] = 10 * (i+1);
	  //cout<<i<<"\t"<<index1[i]<<"\t"<<peakX1[i]<<endl;
	  //cout<<i<<"\t"<<peakX1[i]<<"\t"<<peakX1[i]-peakX1[i-1]<<endl;
	}

      c2->cd(1);
      TGraph *gr1 = new TGraph(ts1->GetNPeaks(),peakX1,time1);
      gr1->SetMarkerStyle(24);
      gr1->Fit("pol1","q");
      fit = gr1->GetFunction("pol1");
      fout<<setw(10)<<fit->GetParameter(0)<<"\t"<<setw(10)<<fit->GetParameter(1)<<"\t";
      gr1->SetTitle(Form("%s Raw_up TDC vs Time ID%d ",detName.Data(),id));
      gr1->GetXaxis()->SetTitle("channel (ch)");
      gr1->GetYaxis()->SetTitle("time (nsec)");
      gr1->Draw("AP");


      c3->cd(1);
      for(Int_t i=0;i<ts1->GetNPeaks();i++)
	res1[i] = time1[i] - fit->Eval(peakX1[i]);
      //cout<<peakX1[ts1->GetNPeaks()-1]<<" "<<time1[ts1->GetNPeaks()-1]<<" "<<fit->Eval(peakX1[ts1->GetNPeaks()-1])<<endl;

      TGraph *gr11 = new TGraph(ts1->GetNPeaks(),peakX1,res1);
      gr11->SetTitle(Form("%s Raw_up Time residual",detName.Data()));
      gr11->GetXaxis()->SetTitle("channel (ch)");
      gr11->GetYaxis()->SetTitle("res (nsec)");
      if(id==1)
	gr11->Draw("AL");
      else
	gr11->Draw("SAMEL");


      
      c1->cd(2);
      tree->Draw("hodTDRaw>>ht2(4950,51.5,5001.5)",Form("hodID==%d",id));
      ht2 = (TH1D*)gDirectory->Get("ht2");
      ht2->SetNameTitle(Form("HOD Raw_down TDC ID%d",id),Form("HOD Raw_down TDC ID%d",id));
      ts2->Search(ht2,2,"",0.1);
      cout<<"ID"<<id<<"(d): "<<ts2->GetNPeaks()<<endl;
      Double_t *peakX2 = new Double_t[ts2->GetNPeaks()];
      Int_t *index2 = new Int_t[ts2->GetNPeaks()];
      Double_t *ch2 = new Double_t[ts2->GetNPeaks()];
      Double_t *err2 = new Double_t[ts2->GetNPeaks()];
      Double_t *time2 = new Double_t[ts2->GetNPeaks()];//psec
      Double_t *res2 = new Double_t[ts2->GetNPeaks()];
      TMath::Sort(ts2->GetNPeaks(),ts2->GetPositionX(),index2,false);
      for(Int_t i=0;i<ts2->GetNPeaks();i++)
	{
	  peakX2[i] = (ts2->GetPositionX())[index2[i]];
	  time2[i] = 10 * (i+1);
	  //cout<<i<<"\t"<<index2[i]<<"\t"<<peakX2[i]<<endl;
	  //cout<<i<<"\t"<<peakX2[i]<<"\t"<<peakX2[i]-peakX2[i-1]<<endl;
	}

      c2->cd(2);
      TGraph *gr2 = new TGraph(ts2->GetNPeaks(),peakX2,time2);
      gr2->SetMarkerStyle(24);
      gr2->Fit("pol1","q");
      fit = gr2->GetFunction("pol1");
      fout<<setw(10)<<fit->GetParameter(0)<<"\t"<<setw(10)<<fit->GetParameter(1)<<endl;
      gr2->SetTitle(Form("%s Raw_down TDC vs Time ID%d ",detName.Data(),id));
      gr2->GetXaxis()->SetTitle("channel (ch)");
      gr2->GetYaxis()->SetTitle("time (nsec)");
      gr2->Draw("AP"); 

      
      c3->cd(2);
      for(Int_t i=0;i<ts2->GetNPeaks();i++)
	res2[i] = time2[i] - fit->Eval(peakX2[i]);
      //cout<<peakX2[ts2->GetNPeaks()-1]<<" "<<time2[ts2->GetNPeaks()-1]<<" "<<fit->Eval(peakX2[ts2->GetNPeaks()-1])<<endl;

      TGraph *gr22 = new TGraph(ts2->GetNPeaks(),peakX2,res2);
      gr22->SetTitle(Form("%s Raw_down Time residual",detName.Data()));
      gr22->GetXaxis()->SetTitle("channel (ch)");
      gr22->GetYaxis()->SetTitle("res (nsec)");
      if(id==1)
	gr22->Draw("AL");
      else
	gr22->Draw("SAMEL");
  

      
      
      c1->Update();
      c2->Update();
      c3->Update();
      //getchar();

      c1->Print("./fig/RawTDC_HOD.pdf");
      c2->Print("./fig/tcal_HOD.pdf");
      c3->Print("./fig/tcal_HOD_res.pdf");

    }
  c1->Print("./fig/RawTDC_HOD.pdf]");
  c2->Print("./fig/tcal_HOD.pdf]");
  c3->Print("./fig/tcal_HOD_res.pdf]");

  fout.close();

  
  /* 
   
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
  */
  //delete [] peakX1,index1,ch1,err1;
  //delete [] peakX2,index2,ch2,err2;

}
