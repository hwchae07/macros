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

TH2 *makeHODQ(Int_t id, Double_t range=10);

void drawHODQ()
{
  TFile *file = new TFile("./hist/empty_run_hodq.root","r");
  TH2 *h1;
  TH1 *hp1,*hp2;
  TF1 *fit1, *fit2;

  //  ofstream fout("./dat/brho_scan_hodq.dat");

  Double_t ratio =1.;

  // fout<<ratio<<std::endl;
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  for(Int_t id=1;id<24;id++)
    {
      file->GetObject(Form("hodq%02d%02d",id,id+1),h1);
      hp1 = h1->ProjectionY(Form("hodq%02d%02d_%02d",id,id+1,id),id,id);
      hp1->SetTitle(Form("HODQ%02d%02d_%02d",id,id+1,id));
      hp2 = h1->ProjectionY(Form("hodq%02d%02d_%02d",id,id+1,id+1),id+1,id+1);
      hp2->SetTitle(Form("HODQ%02d%02d_%02d",id,id+1,id+1));

      cout<<hp1->GetEntries()<<"\t"<<hp2->GetEntries()<<endl;
      /*
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
      */
    }
  //fout.close();
}


void saveHODQ()
{
  TH2 *h1;
  
  for(Int_t id=1;id<24;id++)
    {
      h1 = makeHODQ(id);
      TFile *fhist = new TFile("./hist/empty_run_hodq.root","update");
      h1->Write();
      fhist->Close();
    }

}


TH2 *makeHODQ(Int_t id,Double_t range)
{
  TFile *file3 = new TFile("./cut/pid_empty_run.root","r");
  TCutG *b34na = (TCutG*)file3->Get("b34na");
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

  for(Int_t runNum=295;runNum<=301;runNum++)
    {
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
      //Beam PID//                                                                                                                                              
      Double_t tof713_offset = 213.73055 +352.345;
      Double_t fl7_13 = (3957.1684+3957.2192) / 2.;    //cm                                                                                                     
      Double_t brho0 = 7.8412;    //Tm                                                                                                                          
      tree->SetAlias("tofc713",Form("TOF7_13+%lf",tof713_offset));
      tree->SetAlias("betaF7",Form("%lf/tofc713/29.9792458",fl7_13));
      tree->SetAlias("betaF5","betaF7*1.005");
      tree->SetAlias("betaF13","betaF7*0.9");
      tree->SetAlias("dEFac_beam","TMath::Log(4866 * betaF13**2) - TMath::Log(1-betaF13**2) - betaF13**2");
      tree->SetAlias("beamZ","1+10*TMath::Sqrt(icbEloss/dEFac_beam)*betaF13");
      tree->SetAlias("beamMomU"," 931.494*betaF5/TMath::Sqrt(1-betaF5**2)"); // mom/c = gamma*mass*beta [MeV/u/c]
      tree->SetAlias("beamAoZ",Form("(1+f5X/3300)*%lf * 299.792458 / beamMomU ",brho0));
      //Beam PID//              

      Double_t fdc2Z = 4164.51;    //3726.51 + 438
      Double_t hodZ = 4999.17;    //5004.17-5
      //hodoscope coordinate, 
      //tree->SetAlias("hodX",Form("fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
      tree->SetAlias("hodX",Form("(fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A))*0.9834-5.589 ",hodZ,fdc2Z));

      tree->SetAlias("Qfactor","sqrt(1+tan(fdc2A)**2+tan(fdc2B)**2 )");
      tree->SetAlias("hodQAcor","sqrt(hodQU*hodQD)/Qfactor");
      
      //tree->Draw("hodQAcor:hodID>>hist1(24,0.5,24.5,200,300,600)",Form("b34na&&hodNum==1&&(hodID==%d||hodID==%d)",id,id+1)&&cut_hod,"goff");
      tree->Draw("hodQAcor:hodID>>hist1(24,0.5,24.5,200,300,600)",Form("hodNum==1&&(hodID==%d||hodID==%d)",id,id+1)&&cut_hod,"goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      
    }

  hsum1->SetNameTitle(Form("hodq%02d%02d",id,id+1),Form("hodq%02d%02d",id,id+1));


  return hsum1;

}


