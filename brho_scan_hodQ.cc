#include "TFile.h"
#include "TCutG.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TTree.h"
#include <iostream>
#include <fstream>

TH2 *makeHODQ(Int_t id,Double_t range=10);
TH2 *make2dHODx();
TH2 **make2dHODQx();
Double_t* boundaryHODx();

void drawHODQ()
{
  TFile *file = new TFile("./hist/brho_scan_hodq.root","r");
  TH2 *h1;
  TH1 *hp1,*hp2;
  TF1 *fit1, *fit2;

  ofstream fout("./dat/brho_scan_hodq.dat");

  Double_t ratio =1.;

  fout<<ratio<<std::endl;
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  for(Int_t id=1;id<24;id++)
    {
      file->GetObject(Form("hodq%02d%02d",id,id+1),h1);
      hp1 = h1->ProjectionY(Form("hodq%02d%02d_%02d",id,id+1,id),id,id);
      hp1->SetTitle(Form("HODQ%02d%02d_%02d",id,id+1,id));
      hp2 = h1->ProjectionY(Form("hodq%02d%02d_%02d",id,id+1,id+1),id+1,id+1);
      hp2->SetTitle(Form("HODQ%02d%02d_%02d",id,id+1,id+1));

      c1->cd(1);
      hp1->Draw();
      hp1->Fit("gaus","q");
      fit1 = hp1->GetFunction("gaus");
      hp1->Fit("gaus","rq","",fit1->GetParameter(1)-2*fit1->GetParameter(2),fit1->GetParameter(1)+2*fit1->GetParameter(2));
      fit1 = hp1->GetFunction("gaus");

      c1->cd(2);
      hp2->Draw();
      hp2->Fit("gaus","q");
      fit2 = hp2->GetFunction("gaus");
      hp2->Fit("gaus","rq","",fit2->GetParameter(1)-2*fit2->GetParameter(2),fit2->GetParameter(1)+2*fit2->GetParameter(2));
      fit2 = hp2->GetFunction("gaus");

      ratio *= fit1->GetParameter(1)/fit2->GetParameter(1);
      fout<<ratio<<std::endl;
      std::cout<<fit1->GetParameter(1)/fit2->GetParameter(1)<<std::endl;
      std::cout<<ratio<<std::endl;
      

      
      c1->Update();
      getchar();
    }
  fout.close();
}


void saveHODQ()
{
  TH2 *h1;

  for(Int_t id=1;id<24;id++)
    {
      h1 = makeHODQ(id);
      TFile *fhist = new TFile("./hist/brho_scan_hodq.root","update");
      h1->Write();
      fhist->Close();
    }

}


TH2 *makeHODQ(Int_t id,Double_t range)
{
  TFile *file3 = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne = (TCutG*)file3->Get("b22ne");
  TCut cut_central = "abs(f5X)<2";

  TH2 *hist1;
  TH2 *hsum1 = new TH2D("hsum1","hsum1",24,0.5,24.5,200,300,600);

  ifstream fin("./dat/hodx_boundary.dat");
  Double_t hodx_boundary[24];
  for(Int_t i=1;i<=23;i++)
    fin>>hodx_boundary[i];
  fin.close();
  
  
  TCut cut_hod = Form("abs(hodX-(%lf))<%lf",hodx_boundary[id],range);
  cout<<hodx_boundary[id]<<endl;

  for(Int_t runNum=133;runNum<=138;runNum++)
    {
      //Int_t runNum = 133;
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
      tree->AddFriend("SAMURAI",Form("./root/%s.SAMURAI",filename.Data()));
      tree->GetFriend("SAMURAI")->BuildIndex("EventNum");
  
      Double_t fdc2Z = 4164.51;    //3726.51 + 438
      Double_t hodZ = 4999.17;    //5004.17-5
      //hodoscope coordinate, 
      //tree->SetAlias("hodX",Form("fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
      tree->SetAlias("hodX",Form("(fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A))*0.9834-5.589 ",hodZ,fdc2Z));

      tree->SetAlias("Qfactor","sqrt(1+tan(fdc2A)**2+tan(fdc2B)**2 )");
      tree->SetAlias("hodQAcor","sqrt(hodQU*hodQD)/Qfactor");
      
      tree->Draw("hodQAcor:hodID>>hist1(24,0.5,24.5,200,300,600)",Form("b22ne&&hodNum==1&&(hodID==%d||hodID==%d)",id,id+1)&&cut_hod,"goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      
    }

  hsum1->SetNameTitle(Form("hodq%02d%02d",id,id+1),Form("hodq%02d%02d",id,id+1));


  return hsum1;

}


void test()
{
  Double_t *x = boundaryHODx();
  ofstream fout("./dat/hodx_boundary.dat");

  for(Int_t i=0;i<23;i++)
    fout<<x[i]<<endl;

  fout.close();
}


Double_t* boundaryHODx()
{
  //Double_t *boundary = new Double_t[24];
  
  TH2 *h1 = make2dHODx();
  TH1 *h2;

  TCanvas *c1 = new TCanvas("c1","c1",500,500);

  Double_t hodx_cal[24],hodx_real[24];
  Double_t hodx_min[24],hodx_max[24];
  Double_t *boundary = new Double_t[24];
  
  for(Int_t id=1;id<=24;id++)
    {

      h2 = h1->ProjectionY(Form("HOD%02d",id),id,id);
      //h2->Draw("goff");

      hodx_real[id-1] = 1250.-100.*id;
      hodx_cal[id-1] = h2->GetMean();
      hodx_min[id-1] = h2->GetBinCenter(h2->FindFirstBinAbove(100));
      hodx_max[id-1] = h2->GetBinCenter(h2->FindLastBinAbove(100));

      boundary[id-1] = h2->GetBinCenter(h2->FindFirstBinAbove(100));
      //std::cout<< "ID"<< id <<" "<<1250.-100.*id<<" "<<hodx_min[id-1]<<" "<<h2->GetMean() <<" "<<hodx_max[id-1]<<std::endl;
    }

  return boundary;
  
}



void make1dHODx()
{
  TH2 *h1 = make2dHODx();
  TH1 *h2;

  TCanvas *c1 = new TCanvas("c1","c1",500,500);

  Double_t hodx_cal[24],hodx_real[24];
  Double_t hodx_min[24],hodx_max[24];
  
  for(Int_t id=1;id<=24;id++)
    {

      h2 = h1->ProjectionY(Form("HOD%02d",id),id,id);
      //h2->Draw("goff");

      hodx_real[id-1] = 1250.-100.*id;
      hodx_cal[id-1] = h2->GetMean();
      hodx_min[id-1] = h2->GetBinCenter(h2->FindFirstBinAbove(100));
      hodx_max[id-1] = h2->GetBinCenter(h2->FindLastBinAbove(100));
      std::cout<< "ID"<< id <<" "<<1250.-100.*id<<" "<<hodx_min[id-1]<<" "<<h2->GetMean() <<" "<<hodx_max[id-1]<<std::endl;
      //c1->Update();
      //getchar();
    }
  /*
  TGraph *gr1 = new TGraph(24,hodx_cal,hodx_real);
  gr1->Draw("AP");
  gr1->SetMarkerStyle(20);
  gr1->Fit("pol1");
  */
}

void save2dHODQx()
{
  TH2 **h1 = make2dHODQx();
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  h1[0]->Draw("colz");
  c1->cd(2);
  h1[1]->Draw("colz");
}


TH2** make2dHODQx()
{
  TFile *file3 = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne = (TCutG*)file3->Get("b22ne");
  TCut cut_central = "abs(f5X)<2";

  TH2 **hists = new TH2*[2];
  TH2 *hist1;  
  TH2 *hsum1 = new TH2D("hsum1","hsum1",240,-1200,1200,300,300,600);
  TH2 *hsum2 = new TH2D("hsum2","hsum2",240,-1200,1200,300,300,600);

  for(Int_t runNum=133;runNum<=138;runNum++)
    {
      TString filename = Form("run%04d.root",runNum);
  
      TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("EventNum");
      tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
      tree->GetFriend("FDC")->BuildIndex("EventNum");
      tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
      tree->GetFriend("HOD")->BuildIndex("EventNum");
  
      Double_t fdc2Z = 4164.51;    //3726.51 + 438
      Double_t hodZ = 4999.17;    //5004.17-5
      //hodoscope coordinate, 
      //tree->SetAlias("hodX",Form("fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
      tree->SetAlias("hodX",Form("(fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A))*0.9834-5.589 ",hodZ,fdc2Z));

      tree->SetAlias("Qfactor","sqrt(1+tan(fdc2A)**2+tan(fdc2B)**2 )");
      tree->SetAlias("hodQAcor","sqrt(hodQU*hodQD)/Qfactor");

      tree->SetAlias("hodQcor","hodQ/Qfactor");

      
      tree->Draw("hodQAcor:hodX>>hist1(240,-1200,1200,300,300,600)","b22ne&&hodNum==1","goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);

            
      tree->Draw("hodQcor:hodX>>hist1(240,-1200,1200,300,300,600)","b22ne&&hodNum==1","goff");
      gDirectory->GetObject("hist1",hist1);
      hsum2->Add(hist1);
      
      
    }
  hists[0] = hsum1;
  hists[1] = hsum2;
  return hists;
  
}


TH2 *make2dHODx()
{
  TFile *file3 = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne = (TCutG*)file3->Get("b22ne");
  TCut cut_central = "abs(f5X)<2";

  TH2 *hist1;
  TH2 *hsum1 = new TH2D("hsum1","hsum1",24,0.5,24.5,5000,-1300,1300);

  for(Int_t runNum=133;runNum<=138;runNum++)
    {
      TString filename = Form("run%04d.root",runNum);
  
      TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("EventNum");
      tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
      tree->GetFriend("FDC")->BuildIndex("EventNum");
      tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
      tree->GetFriend("HOD")->BuildIndex("EventNum");
  
      Double_t fdc2Z = 4164.51;    //3726.51 + 438
      Double_t hodZ = 4999.17;    //5004.17-5
      //hodoscope coordinate, 
      //tree->SetAlias("hodX",Form("fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
      tree->SetAlias("hodX",Form("(fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A))*0.9834-5.589 ",hodZ,fdc2Z));
      
      tree->Draw("hodX:hodID>>hist1(24,0.5,24.5,200,-1300,1300)","b22ne&&hodNum==1","goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      
    }

  return hsum1;

}
