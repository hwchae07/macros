{
  Int_t runNum=297;
  
  TFile *file = new TFile(Form("./root/run%04d.root.BeamPla",runNum));
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("EventNum");

  tree->AddFriend("BDC",Form("./root/run%04d.root.BDC",runNum));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  tree->AddFriend("HOD",Form("./root/run%04d.root.HOD",runNum));
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  tree->AddFriend("FDC",Form("./root/run%04d.root.FDC",runNum));
  tree->GetFriend("FDC")->BuildIndex("EventNum");


  
  
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
  Double_t tgtZ = -4978.89+488;
	     
  tree->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
  tree->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
  TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";  
  //Target size//

  Double_t fdc2Z = 4164.51;    //3726.51 + 438
  Double_t hodZ = 4999.17;    //5004.17-5
  tree->SetAlias("hodX",Form("(fdc2X + 715.01 - 813.88  + (%lf-%lf)*tan(fdc2A))*0.9834-5.589 ",hodZ,fdc2Z));

  tree->SetAlias("Qfactor","sqrt(1+tan(fdc2A)**2+tan(fdc2B)**2 )");
  tree->SetAlias("hodQAcor","sqrt(hodQU*hodQD)/Qfactor");
  tree->SetAlias("hodQcor","hodQ/Qfactor");

  /*
  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  tree->Draw("hodQAcor:hodX>>h1(500,100,300,500,300,500)","","colz");
  c1->cd(2);
  tree->Draw("hodQcor:hodX>>h2(500,100,300,500,300,500)","","colz");
  */

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  tree->Draw("hodQcor:hodX>>h1(2000,-1000,1000,500,200,600)","","colz");     
  
}
