{

  TH1 *hodq1 = new TH1D("hodq1","hodq;hodQ",1000,200,550);
  TH1 *hodqsub1 = new TH1D("hodqsub1","hodq;hodQ",1000,200,550);

  TH2 *hsum = new TH2D("hsum","hsum;FDC2X;hodT",200,-1200,500,200,80,90);
  TH2 *hist = new TH2D("hist","hist",200,-1200,500,200,80,90);
  TH2 *hsumQ1 = new TH2D("hsumQ1","hsumQ;hodT;hodQ",200,30,65,200,200,550);
  TH2 *histQ1 = new TH2D("histQ1","histQ",200,30,65,200,200,550);

  TH2 *hsumPID = new TH2D("hsumPID","Z vs A/Z;A/Z;Z",200,2.0,3.5,200,3.5,5);
  TH2 *histPID = new TH2D("histPID","Z vs A/Z;A/Z;Z",200,2.0,3.5,200,3.5,5);
  
  TCanvas *c2 = new TCanvas("c2","c2",800,800);

  TFile *fcut = new TFile("./cut/pid_32ne_physics.root","r");
  TCutG *f32ne;
  fcut->GetObject("f32ne",f32ne);
  f32ne->SetVarX("fragAoZ");
  f32ne->SetVarY("fragZ");
  //  c2->Divide(2,1);

  Int_t runNum = 275;
  //for(Int_t runNum=275 ; runNum<=294 ; runNum++)
  //{
  
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
  tree->AddFriend("Coin",Form("./root/%s.Coin",filename.Data()));
  tree->GetFriend("Coin")->BuildIndex("EventNum");
  
  
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

  //HOD time offset//
  TString tof_hod = "";

  tof_hod += "hodTA-T13";
  tree->SetAlias("tof_hod",tof_hod);
  //HOD time offset//

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
 
  tree->SetAlias("hodzz",Form("%lf*cos(%lf) - (hodX)*sin(%lf)",hodZ,alpha,alpha));
  tree->SetAlias("hodxx",Form("%lf*sin(%lf) + (hodX)*cos(%lf)",hodZ,alpha,alpha));
  tree->SetAlias("theta",Form("abs(%lf + fdc1A - fdc2A)",alpha));
      
  tree->SetAlias("FL1",Form("sqrt( (interX-tgtX)**2 + (interZ-%lf)**2 )",tgtZ));
  tree->SetAlias("FL2","sqrt( (interX-hodxx)**2 + (interZ-hodzz)**2 )");
  tree->SetAlias("FL","FL1 + FL2 + 1000*brho/2.9*theta - 2*1000*brho/2.9*tan(theta/2)");


  tree->SetAlias("beta","FL/tof_hod/299.792458");


  tree->SetAlias("dEFac","TMath::Log(15795.98 * beta * beta) - TMath::Log(1-beta*beta) - beta*beta");
  tree->SetAlias("fragZ","TMath::Sqrt(hodQ/TMath::Sqrt( (1+tan(fdc2A)**2 + tan(fdc2B)**2) )/dEFac)*beta");

  tree->SetAlias("fragMomUfromHOD"," 931.494*beta/TMath::Sqrt(1-beta*beta)"); // mom/c = gamma*mass*beta [MeV/u/c]
  tree->SetAlias("fragAoZ","brho * 299.792458 / fragMomUfromHOD ");


  //coincidence with neutron//
  tree->SetAlias("trig_neut","coinTrigger==3||coinTrigger==5||coinTrigger==7");
  //coincidence with neutron//
      
  
  
  
  /*
    c2->cd(1);
    tree->Draw("hodQA>>hodqsub1(1000,200,550)",cut_32ne&&cut_reaction);
    gDirectory->GetObject("hodqsub1",hodqsub1);     
    hodq1->Add(hodqsub1);
    gPad->SetLogy();
    hodq1->Draw();

    c2->cd(2);
    tree->Draw("hodQA:tof_hod>>histQ1(200,30,65,200,200,550)",cut_32ne&&cut_reaction,"col");
    gDirectory->GetObject("histQ1",histQ1);
    hsumQ1->Add(histQ1);
    hsumQ1->Draw("col");
  */

  tree->Draw("fragZ:fragAoZ>>histPID(200,2.0,3.5,200,3.5,5)",cut_34na&&cut_reaction&&"hodNum==1","colz");
  gDirectory->GetObject("histPID",histPID);
  hsumPID->Add(histPID);
  hsumPID->Draw("colz");
  c2->Update();
  //}

}
