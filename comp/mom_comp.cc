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

  Int_t runNum = 300;
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
  tree->AddFriend("FragPID",Form("./root/%s.FragPID",filename.Data()));
  tree->GetFriend("FragPID")->BuildIndex("EventNum");

  //physics run
  //TCut f34na = "abs(fragAoZ-3.17)<0.05&&abs(fragZ-4.2)<0.15";
  //physics run
  //empty run
  TCut f34na = "abs(fragAoZ-3.24)<0.04&&abs(fragZ-4.3)<0.1";
  //emptr run
  //tree->Draw("fragZ:fragAoZ>>h2(400,2.0,3.5,400,3.5,5.0)","b34na"&&f34na,"colz");

  TCanvas *c1 = new TCanvas("c1","c1",1200,400);
  c1->Divide(3,1);
  c1->cd(1);
  tree->Draw("beamAoZ:fdc2X>>h1(200,150,450,200,2.7,3.3)","b34na","colz");
  c1->cd(2);
  tree->Draw("brho*299.792458/beamMomU:fdc2X>>h2(200,150,450,200,2.7,3.3)","b34na","colz");
  c1->cd(3);
  tree->Draw("brho*299.792458/(941.494/sqrt(1-betaFrag*betaFrag)*betaFrag):fdc2X>>h3(200,150,450,200,2.7,3.3)","b34na","colz");
  
}
