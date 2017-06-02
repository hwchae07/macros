{
  gROOT->Reset();

  Int_t runNum = 133;
  
  TString filename = Form("run%04d.root",runNum);
  
  TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("EventNum");
  tree->AddFriend("BDC",Form("./root/%s.BDC",filename.Data()));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  //tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
  //tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  //tree->AddFriend("Brho",Form("./root/%s.Brho",filename.Data()));
  //tree->GetFriend("Brho")->BuildIndex("EventNum");
  //tree->AddFriend("SAMURAI",Form("./root/%s.SAMURAI",filename.Data()));
  //tree->GetFriend("SAMURAI")->BuildIndex("EventNum");
  //tree->AddFriend("Coin",Form("./root/%s.Coin",filename.Data()));
  //tree->GetFriend("Coin")->BuildIndex("EventNum");

  TFile *file3 = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne = (TCutG*)file3->Get("b22ne");
  
  TCut cut_central = "abs(f5X)<5";
  
  //TOF F7~F13 offset setting//
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  tree->Draw("TOF7_13>>h1(100,-368,-366)",cut_central&&"b22ne");
  
  TH1 *h1;
  gDirectory->GetObject("h1",h1);
  TF1 *fit1 = new TF1("fit1","gaus",-368,-366);
  h1->Fit(fit1,"rq","",-367,-366.5);


  Double_t offset_pla = 197.78065 - fit1->GetParameter(1);
  
  tree->SetAlias("tofc713",Form("TOF7_13+%lf",offset_pla));

  Double_t par[4] = {1.11972578125000000000E+04,
		     -1.41039001464843750000E+02,
		     6.17486000061035156250E-01,
		     -9.23150219023227691650E-04};

  //TOF F7~F13 to Energy//
  TString energy = "";
  energy += Form("%lf",par[0]);
  for(Int_t i=1;i<4;i++)
    energy += Form("+ (%lf) * TMath::Power(tofc713,%d)",par[i],i);

  Double_t fl7_13 = (3957.1684+3957.2192) / 2. ;    //cm

  tree->SetAlias("energy",energy);

  
}
