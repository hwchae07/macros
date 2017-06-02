{
  //PID cut//
  TFile *pidcut = new TFile("./cut/pid_32ne_physics.root","r");
  TCutG *b34na = (TCutG*)pidcut->Get("b34na");
  TCutG *f32ne = (TCutG*)pidcut->Get("f32ne");
  f32ne->SetVarX("fragAoZ");
  f32ne->SetVarY("fragZ");
  //PID cut//

  Int_t runNum=275;

  TFile *file = new TFile(Form("./root/run%04d.root.BeamPla",runNum),"r");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("EventNum");
  tree->AddFriend("SAMURAI",Form("./root/run%04d.root.SAMURAI",runNum));
  tree->GetFriend("SAMURAI")->BuildIndex("EventNum");
  tree->AddFriend("FDC",Form("./root/run%04d.root.FDC",runNum));
  tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->AddFriend("BDC",Form("./root/run%04d.root.BDC",runNum));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  tree->AddFriend("HOD",Form("./root/run%04d.root.HOD",runNum));
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  tree->AddFriend("Coin",Form("./root/run%04d.root.Coin",runNum));
  tree->GetFriend("Coin")->BuildIndex("EventNum");

  //Beam PID//
  Double_t tof713_offset = 213.73055 +352.345;
  Double_t fl7_13 = (3957.1684+3957.2192) / 2.;    //cm
  Double_t brho0 = 7.8412;    //Tm
  tree->SetAlias("tofc713",Form("TOF7_13+%lf",tof713_offset));
  tree->SetAlias("betaF7",Form("%lf/tofc713/29.9792458",fl7_13));
  tree->SetAlias("betaF5","betaF7*1.005");
  tree->SetAlias("betaF13","betaF7*0.9");
  tree->SetAlias("dEFac_beam","log(4866 * betaF13 * betaF13) - log(1-betaF13*betaF13) - betaF13*betaF13");
  tree->SetAlias("beamZ","1+10*TMath::Sqrt(icbEloss/dEFac_beam)*betaF13");
  tree->SetAlias("beamMomU"," 931.494*betaF5/TMath::Sqrt(1-betaF5*betaF5)"); // mom/c = gamma*mass*beta [MeV/u/c]
  tree->SetAlias("beamAoZ",Form("(1+f5X/3300)*%lf * 299.792458 / beamMomU ",brho0));
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
  //ifstream file1("./dat/hod_offset.dat");
  ifstream file1("./dat/hod_offset_brho_scan.dat");
  Int_t dummy[25];
  Double_t offset[25];
  for(Int_t i=0;i<=24;i++)
    file1>>dummy[i]>>offset[i];
  TString tof_hod = "";
  tof_hod += "hodTA-T13";
  tof_hod += "+";
  tof_hod += offset[0];
  for(Int_t id=1;id<=24;id++)
    {
      tof_hod += "+";
      tof_hod += offset[id];
      tof_hod += Form("*(hodID==%d)",id);
    }
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

  //hod position//
  Double_t hodZ = 4999.17;
  tree->SetAlias("hodX",Form("fdc2X + 715.01 + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
  tree->SetAlias("hodzz",Form("%lf*cos(%lf) - (hodX)*sin(%lf)",hodZ,alpha,alpha));
  tree->SetAlias("hodxx",Form("%lf*sin(%lf) + (hodX)*cos(%lf)",hodZ,alpha,alpha));
  tree->SetAlias("theta",Form("abs(%lf + fdc1A - fdc2A)",alpha));
  //hod position//
  
  //fragment PID//
  tree->SetAlias("FL1",Form("sqrt( (interX-tgtX)**2 + (interZ-%lf)**2 )",tgtZ));
  tree->SetAlias("FL2","sqrt( (interX-hodxx)**2 + (interZ-hodzz)**2 )");
  tree->SetAlias("FL","FL1 + FL2 + 1000*brho/2.9*theta - 2*1000*brho/2.9*tan(theta/2)");

  tree->SetAlias("beta","FL/tof_hod/299.792458");
  tree->SetAlias("dEFac_frag","log(15795.98 * beta * beta) - log(1-beta*beta) - beta*beta");
  tree->SetAlias("fragZ","sqrt(hodQ/sqrt(1 + tan(fdc2A)**2 + tan(fdc2B)**2)/dEFac_frag)*beta");
  tree->SetAlias("fragMomUfromHOD"," 931.494*beta/TMath::Sqrt(1-beta*beta)"); // mom/c = gamma*mass*beta [MeV/u/c]
  tree->SetAlias("fragAoZ","brho * 299.792458 / fragMomUfromHOD ");
  //fragment PID//
  
  TCut fragbox = "abs(fragZ-4)<1&&abs(fragAoZ-2.7)<0.7";
  
  
}
