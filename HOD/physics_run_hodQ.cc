#include "TFile.h"
#include "TCutG.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TTree.h"
#include "TMath.h"
#include <iostream>
#include <fstream>

TH2 *makeHODQ(Int_t id, Double_t range=10);

void drawHODQ()
{
  TFile *file = new TFile("./hist/physics_run_hodq.root","r");
  TH2 *h1;
  TH1 *hp1,*hp2;
  TF1 *fit1, *fit2;

  ofstream fout("./dat/physics_run_hodq.dat");

  Double_t ratio =1.;

  fout<<ratio<<std::endl;
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  c1->Print("./fig/physics_run_hodQ.pdf[");
  for(Int_t id=8;id<24;id++)
    {
      file->GetObject(Form("hodq%02d%02d",id,id+1),h1);
      hp1 = h1->ProjectionY(Form("hodq%02d%02d_%02d",id,id+1,id),id,id);
      hp1->SetTitle(Form("HODQ%02d%02d_%02d",id,id+1,id));
      hp1->Rebin(5);
      hp2 = h1->ProjectionY(Form("hodq%02d%02d_%02d",id,id+1,id+1),id+1,id+1);
      hp2->SetTitle(Form("HODQ%02d%02d_%02d",id,id+1,id+1));
      hp2->Rebin(5);
      
      //cout<<hp1->GetEntries()<<"\t"<<hp2->GetEntries()<<endl;

      TSpectrum *ts1 = new TSpectrum();
      TSpectrum *ts2 = new TSpectrum();


      
      c1->cd(1);
      hp1->Draw();
      ts1->Search(hp1,2,"",0.5);
      Int_t nPeaks1 = ts1->GetNPeaks();
      Double_t *x1 = new Double_t[nPeaks1];
      Int_t *index1 = new Int_t[nPeaks1];
      x1 = ts1->GetPositionX();
      TMath::Sort(nPeaks1,x1,index1);
      for(Int_t i=0;i<nPeaks1;i++)
	x1[i] = (ts1->GetPositionX())[index1[i]];
      Double_t center1 = x1[0];
      Double_t range1 = 20;
      hp1->Fit("gaus","rq","",center1-range1,center1+range1);
      fit1 = hp1->GetFunction("gaus");
      //cout<<center1<<","<<fit1->GetParameter(1)<<endl;
      hp1->Fit("gaus","rq","",fit1->GetParameter(1)-2*fit1->GetParameter(2),fit1->GetParameter(1)+2*fit1->GetParameter(2));
      fit1 = hp1->GetFunction("gaus");

      c1->cd(2);
      ts2->Search(hp2,2,"",0.5);
      Int_t nPeaks2 = ts2->GetNPeaks();
      Double_t *x2 = new Double_t[nPeaks2];
      Int_t *index2 = new Int_t[nPeaks2];
      x2 = ts2->GetPositionX();
      TMath::Sort(nPeaks2,x2,index2);
      for(Int_t i=0;i<nPeaks2;i++)
	x2[i] = (ts2->GetPositionX())[index2[i]];
      hp2->Draw();
      Double_t center2 = x2[0];
      Double_t range2 = 20;
      hp2->Fit("gaus","rq","",center2-range2,center2+range2);
      fit2 = hp2->GetFunction("gaus");
      hp2->Fit("gaus","rq","",fit2->GetParameter(1)-2*fit2->GetParameter(2),fit2->GetParameter(1)+2*fit2->GetParameter(2));
      fit2 = hp2->GetFunction("gaus");

      /*
	cout<<x1[0]<<"\t"<<x2[0]<<endl;
	cout<<ts1->GetNPeaks()<<"\t"<<ts2->GetNPeaks()<<endl;
      */
      
      ratio *= fit1->GetParameter(1)/fit2->GetParameter(1);
      fout<<ratio<<std::endl;
      std::cout<<fit1->GetParameter(1)/fit2->GetParameter(1)<<std::endl;
      std::cout<<ratio<<std::endl;
     

      
      c1->Update();
      c1->Print("./fig/physics_run_hodQ.pdf");
      getchar();
      
    }
    c1->Print("./fig/physics_run_hodQ.pdf]");
  fout.close();
}


void saveHODQ()
{
  TH2 *h1;
  
  for(Int_t id=1;id<24;id++)
    {
      h1 = makeHODQ(id);
      TFile *fhist = new TFile("./hist/physics_run_hodq.root","update");
      h1->Write();
      delete fhist;
    }

}


TH2 *makeHODQ(Int_t id,Double_t range)
{
  TFile *file3 = new TFile("./cut/pid_32ne_physics.root","r");
  TCutG *b34na = (TCutG*)file3->Get("b34na");
  TCut cut_na = "abs(elossc-22.5)<2.5";
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

  for(Int_t runNum=275;runNum<=294;runNum++)
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

      tree->SetAlias("tofc","TOF3_13-3*dT5");
      tree->SetAlias("elossc","icbEloss-0.1*tofc-40");

                  
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
      tree->Draw("hodQAcor:hodID>>hist1(24,0.5,24.5,200,300,600)",Form("hodNum==1&&(hodID==%d||hodID==%d)",id,id+1)&&cut_hod&&cut_na,"goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      
      delete file;
    }

  hsum1->SetNameTitle(Form("hodq%02d%02d",id,id+1),Form("hodq%02d%02d",id,id+1));

  //delete file3;
    
  return hsum1;

  
}


void make2dHODQ()
{
  TFile *file3 = new TFile("./cut/pid_32ne_physics.root","r");
  TCutG *b34na = (TCutG*)file3->Get("b34na");
  TCut cut_na = "abs(elossc-22.5)<2.5";
  TCut cut_central = "abs(f5X)<2";

  TH2 *hist1;
  TH2 *hsum1 = new TH2D("hsum1","hsum1",24,0.5,24.5,200,300,600);
  TH2 *hsum2 = new TH2D("hsum2","hsum2",2000,-1000,1000,200,300,600);

  ifstream fin("./dat/hodx_boundary.dat");
  Double_t hodx_boundary[24];
  for(Int_t i=1;i<=23;i++)
    fin>>hodx_boundary[i];
  fin.close();
  
  
  //TCut cut_hod = Form("abs(hodX-(%lf))<%lf",hodx_boundary[id],range);
  //cout<<hodx_boundary[id]<<endl;

  for(Int_t runNum=275;runNum<=294;runNum++)
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

      tree->SetAlias("tofc","TOF3_13-3*dT5");
      tree->SetAlias("elossc","icbEloss-0.1*tofc-40");

                  
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
      tree->SetAlias("hodQcor","hodQ/Qfactor");
      
      //tree->Draw("hodQAcor:hodID>>hist1(24,0.5,24.5,200,300,600)",Form("b34na&&hodNum==1&&(hodID==%d||hodID==%d)",id,id+1)&&cut_hod,"goff");
      tree->Draw("hodQAcor:hodID>>hist1(24,0.5,24.5,200,300,600)",cut_na,"goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);

      //tree->Draw("hodQAcor:hodX>>hist1(24,0.5,24.5,200,300,600)","hodNum==1"&&cut_na,"goff");
      //gDirectory->GetObject("hist1",hist1);
      //hsum1->Add(hist1);

      delete file;
    }

  TCanvas *c1 = new TCanvas("c1","c1",600,400);
  hsum1->SetTitle("HODF24 Q VS ID;ID;Q (arbitrary)");
  hsum1->Draw("colz");

  
}


TH2 *makePIDfrag(Bool_t align)
{

  TFile *fcut = new TFile("./cut/beamCut.root","r");
  TCutG *b32na = (TCutG*)fcut->Get("b32na");
  TCutG *b33na = (TCutG*)fcut->Get("b33na");
  TCutG *b34na = (TCutG*)fcut->Get("b34na");
  TCutG *b30ne = (TCutG*)fcut->Get("b30ne");
  TCutG *b31ne = (TCutG*)fcut->Get("b31ne");
  TCutG *b32ne = (TCutG*)fcut->Get("b32ne");
  TCutG *b29f = (TCutG*)fcut->Get("b29f");
  fcut->Close();

  TH2 *hist1;
  TH2 *hsum1 = new TH2D("hsum1","hsum1",400,2.0,3.5,400,3.5,5.0);
  //  tree->Draw("fragZ:fragAoZ>>hist1(200,2.0,3.5,200,3.5,5.0)","b34na","goff");

  for(Int_t runNum=275;runNum<=294;runNum++)
    {
      TString filename = Form("run%04d.root",runNum);
  
      TFile *file = new TFile(Form("./root/%s.BeamPID",filename.Data()),"r");
      TTree *tree = (TTree*)file->Get("BeamPID");
      tree->BuildIndex("EventNum");
      tree->AddFriend("BDC",Form("./root/%s.BDC",filename.Data()));
      tree->GetFriend("BDC")->BuildIndex("EventNum");
      tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
      tree->GetFriend("FDC")->BuildIndex("EventNum");
      tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
      tree->GetFriend("HOD")->BuildIndex("EventNum");
      tree->AddFriend("SAMURAI",Form("./root/%s.SAMURAI",filename.Data()));
      tree->GetFriend("SAMURAI")->BuildIndex("EventNum");


      //Target size//                                                                                                             
      Double_t bdcWidth = 120.;
      Double_t bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
      Double_t bdc1Tgt = bdc2Tgt + 1000;
      Double_t tgtZ = -4978.89+488;

      tree->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
      tree->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
      TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";
      //Target size//                                                                                                             

      TString tof_hod = "";
      tof_hod += "hodTA-T13";
      tree->SetAlias("tof_hod",tof_hod);

      //FDC z position//                                                                                                          
      Double_t fdc1Z = -2888.82;
      Double_t alpha = -59.92/180.*TMath::Pi();
      Double_t fdc2Z = 4164.51;    //3726.51+438;                                                                                 

      tree->SetAlias("fdc2zz",Form("%lf*cos(%lf) - (fdc2X+715.01)*sin(%lf)",fdc2Z,alpha,alpha));
      tree->SetAlias("fdc2xx",Form("%lf*sin(%lf) + (fdc2X+715.01)*cos(%lf)",fdc2Z,alpha,alpha));
      //FDC z position//                                                                                                          

      //intersection//                                                                                                            
      tree->SetAlias("a1",Form("tan( fdc1A+%lf )",0.));
      tree->SetAlias("b1",Form("fdc1X - %lf*a1",fdc1Z));
      tree->SetAlias("a2",Form("tan(fdc2A+%lf)",alpha));
      tree->SetAlias("b2","fdc2xx - fdc2zz*a2");
      tree->SetAlias("interZ","(b2-b1)/(a1-a2)");
      tree->SetAlias("interX","(a1*b2-b1*a2)/(a1-a2)");
      //intersection//                                                                                                            


      Double_t hodZ = 4999.17;
      tree->SetAlias("hodX",Form("fdc2X + 715.01 + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
      //tree->SetAlias("hodX","2667.14-2667.14/24/2-hodID*2667.14/24");                                                           
      tree->SetAlias("hodzz",Form("%lf*cos(%lf) - (hodX)*sin(%lf)",hodZ,alpha,alpha));
      tree->SetAlias("hodxx",Form("%lf*sin(%lf) + (hodX)*cos(%lf)",hodZ,alpha,alpha));
      tree->SetAlias("theta",Form("abs(%lf + fdc1A - fdc2A)",alpha));

      tree->SetAlias("FL1",Form("sqrt( (interX-tgtX)**2 + (interZ-%lf)**2 )",tgtZ));
      tree->SetAlias("FL2","sqrt( (interX-hodxx)**2 + (interZ-hodzz)**2 )");
      tree->SetAlias("FL","FL1 + FL2 + 1000*brho/2.9*theta - 2*1000*brho/2.9*tan(theta/2)");

      tree->SetAlias("beta","FL/tof_hod/299.792458");



      tree->SetAlias("dEFac","TMath::Log(15795.98 * beta * beta) - TMath::Log(1-beta*beta) - beta*beta");

      if(align==true)
	tree->SetAlias("fragZ","TMath::Sqrt(hodQ/TMath::Sqrt( (1+tan(fdc2A)**2 + tan(fdc2B)**2) )/dEFac)*beta");
      else
	tree->SetAlias("fragZ","TMath::Sqrt(hodQA/TMath::Sqrt( (1+tan(fdc2A)**2 + tan(fdc2B)**2) )/dEFac)*beta");
	

      tree->SetAlias("fragMomUfromHOD"," 931.494*beta/TMath::Sqrt(1-beta*beta)"); // mom/c = gamma*mass*beta [MeV/u/c]            
      tree->SetAlias("fragAoZ","brho * 299.792458 / fragMomUfromHOD ");


      //tree->Draw("fragZ:fragAoZ>>histPID(200,2.0,3.5,200,3.5,5)",cut_34na&&cut_reaction&&"hodNum==1&&trig_neut","colz");
      tree->Draw("fragZ:fragAoZ>>hist1(400,2.0,3.5,400,3.5,5.0)","b34na","goff");
      
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      
      delete file;
    }

  //hsum1->SetNameTitle(Form("hodq%02d%02d",id,id+1),Form("hodq%02d%02d",id,id+1));

  //delete file3;
    
  return hsum1;

  
}


void savePIDfrag()
{
  TH2 *h1;
  TH2 *h2;

  h1 = makePIDfrag(true);
  h2 = makePIDfrag(false);
  h1->SetNameTitle("h1","afterAlign");
  h2->SetNameTitle("h2","beforeAlign");
  TFile *fhist = new TFile("./hist/physics_run_fragPID.root","recreate");
  h1->Write();
  h2->Write();
  delete fhist;


}
