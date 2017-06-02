#include "TFile.h"
#include "TTree.h"

void beamPID()
{
  TFile *file = new TFile("root/run0275.root.BeamPla","read");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->SetAlias("tofc","TOF3_13-3*dT5");
  tree->SetAlias("elossc","icbEloss-0.1*tofc-40");

  TCut cut_na = "abs(elossc-22.5)<2.5";
  TCut cut_ne = "abs(elossc-17.8)<2.2";
  TCut cut_32ne = cut_ne&&"abs(tofc+379)<3";
  TCut cut_34na = cut_na&&"abs(tofc+388)<3";

  
  TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  c1->Divide(2,2);
  c1->cd(1);
  //tree->Draw("icbEloss>>heloss(300,10,40)");
  tree->Draw("elossc>>heloss(300,10,30)","");

  c1->cd(2);
  //tree->Draw("TOF3_13:dT5>>hist1(200,-5,-1.5,200,-420,-370)","","col");
  //tree->Draw("icbEloss:TOF3_13-3*dT5>>hist1(200,-420,-370,200,10,30","","col");
  tree->Draw("elossc:tofc>>hist1(200,-410,-360,200,10,30","","col");

  c1->cd(3);
  tree->Draw("elossc>>hna(300,10,30)",cut_na);

  c1->cd(4);
  tree->Draw("elossc:tofc>>hist2(200,-410,-360,200,10,30",cut_34na,"col");
  

  ofstream out("./dat/beamNumber.dat");
  
  TFile *f1;
  TTree *t1;
  TH2 *hist;
  Int_t id;
  Double_t num_tot,num_34na,num_34na_target,ratio,ratio_target;
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  c2->Divide(2,2);
  c2->Print("./fig/beamPID.pdf[");


  out<<"RunNum\tTotal\t34Na\tratio\t34Na_re\tratio_re"<<endl;
  
		      
  for(Int_t runNum=275;runNum<=294;runNum++)
    {
      id = runNum-274;
      c2->cd((id-1)%4+1);
      f1 = new TFile(Form("root/run%04d.root.BeamPla",runNum),"read");      
      f1->GetObject("BeamPla",t1);
      t1->BuildIndex("EventNum");
      t1->AddFriend("BDC",Form("./root/run%04d.root.BDC",runNum));
      t1->GetFriend("BDC")->BuildIndex("EventNum");
      
      t1->SetAlias("tofc","TOF3_13-3*dT5");
      t1->SetAlias("elossc","icbEloss-0.1*tofc-40");

      //Target size//
      Double_t bdcWidth = 120.;
      Double_t bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
      Double_t bdc1Tgt = bdc2Tgt + 1000;
      t1->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
      t1->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
      TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";
      //Target size//


      //cout<<t1<<endl;
      t1->Draw(Form("elossc:tofc>>h%d(200,-410,-360,200,10,30",runNum),cut_34na,"col");
      gDirectory->GetObject(Form("h%d",runNum),hist);

      num_tot = t1->GetEntries();
      num_34na = t1->GetEntries(cut_34na);
      num_34na_target = t1->GetEntries(cut_34na&&cut_reaction);
      ratio = num_34na / num_tot * 100;
      ratio_target = num_34na_target / num_tot *100;
      out<<"run0"<<runNum<<"\t"<<num_tot<<"\t"<<num_34na<<"\t"<< std::setw(7) << ratio<<"\t"<<num_34na_target<<"\t"<<ratio_target<<endl;
      hist->SetDirectory(0);

      c2->Update();

      if(id%4==0)
	c2->Print("./fig/beamPID.pdf");
      
      f1->Close();
      delete f1;
    }
  c2->Print("./fig/beamPID.pdf]");
  
}
