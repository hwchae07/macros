{
  TFile *file = new TFile("root/run0275.root.BeamPla","read");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("RunNum","EventNum");


  tree->SetAlias("tofc","TOF3_13-3*dT5");
  tree->SetAlias("elossc","icbEloss-0.1*tofc-40");

  TCut box = "abs(tofc+360)<50&&abs(elossc-18)<8";
  
  TCut cut_na = "abs(elossc-22.5)<2.5";
  TCut cut_ne = "abs(elossc-17.8)<2.2";
  TCut cut_32ne = cut_ne&&"abs(tofc+379)<3";
  TCut cut_34na = cut_na&&"abs(tofc+388)<3";

  TCut cut_central = "abs(f5X)<10";
  
  tree->Draw("TOF7_13>>h1(100,-354,-351)",cut_32ne&&cut_central);
  TH1 *h1;
  gDirectory->GetObject("h1",h1);
  TF1 *fit1 = new TF1("fit1","gaus",-354,-351);
  h1->Fit(fit1,"rq","",-353,-351.5);

  cout<<fit1->GetParameter(1)<<endl;
  Double_t offset = 213.73055 - fit1->GetParameter(1);

  tree->SetAlias("tofc713",Form("TOF7_13+%lf",offset));
  //tree->Draw("tofc713>>(100,200,220)");
  
  
  //Parameters for TOF(F7-F13) to Energy(B.Tg)
  Double_t par[4] = { 6.84494091796875000000E+03,
		      -7.73498382568359375000E+01,
		      3.06769594550132751465E-01,
		      -4.17593779275193810463E-04};
  TString energy = "";
  energy += Form("%lf",par[0]);
  for(Int_t i=1;i<4;i++)
    energy += Form("+ (%lf) * TMath::Power(tofc713,%d)",par[i],i);

  Double_t fl7_13 = (3957.1684+3957.2192) / 2. ;    //cm
  Double_t brho0 = 7.8412;    //Tm

  tree->SetAlias("beta",Form("%lf/tofc713/29.9792458",fl7_13));

  tree->SetAlias("betaF5","beta*1.005");
  tree->SetAlias("betaF13","beta*0.9");
  
  tree->SetAlias("dEFac","TMath::Log(4866 * betaF13 * betaF13) - TMath::Log(1-betaF13*betaF13) - betaF13*betaF13");
  tree->SetAlias("beamZ","1+10*TMath::Sqrt(icbEloss/dEFac)*betaF13");

  tree->SetAlias("beamMomU"," 931.494*betaF5/TMath::Sqrt(1-betaF5*betaF5)"); // mom/c = gamma*mass*beta [MeV/u/c]
  tree->SetAlias("beamAoZ",Form("(1+f5X/3300)*%lf * 299.792458 / beamMomU ",brho0));

  tree->Draw("beamZ:beamAoZ>>hi1(400,2.8,3.5,400,8,12)","","colz");
  TH1 *hi1 = (TH1D*)gDirectory->Get("hi1");

  

  TFile *file2 = new TFile("./root/run0275.root.BeamPID","r");
  TTree *tree2 = (TTree*)file2->Get("BeamPID");
  tree2->BuildIndex("RunNum","EventNum");
  tree2->AddFriend("BeamPla","./root/run0275.root.BeamPla");
  tree2->GetFriend("BeamPla")->BuildIndex("RunNum","EventNum");

  tree2->SetAlias("tofc","TOF3_13-3*dT5");
  tree2->SetAlias("elossc","icbEloss-0.1*tofc-40");
  
  tree2->Draw("beamZ:beamAoZ>>hi2(400,2.8,3.5,400,8,12)","","colz");
  TH1 *hi2 = (TH1D*)gDirectory->Get("hi2");
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  hi1->Draw("colz");
  hi1->SetTitle("Z VS A/q;A/q;Z");
  c1->cd(2);
  hi2->Draw("colz");
  hi2->SetTitle("Z VS A/q;A/q;Z");

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  c2->cd(1);
  tree->Draw("beamAoZ>>hi3(400,2.8,3.5)",cut_na);
  hi3->SetTitle("Z VS A/q;A/q;Counts");
  gPad->SetLogy();
  c2->cd(2);
  tree2->Draw("beamAoZ>>hi4(400,2.8,3.5)",cut_na);
  gPad->SetLogy();
  hi4->SetTitle("Z VS A/q;A/q;Counts");
}
