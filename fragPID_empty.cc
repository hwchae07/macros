{

  gStyle->SetOptStat(0);
  
  TFile *fcut = new TFile("./cut/beamCut.root","r");
  TCutG *b32na = (TCutG*)fcut->Get("b32na");
  TCutG *b33na = (TCutG*)fcut->Get("b33na");
  TCutG *b34na = (TCutG*)fcut->Get("b34na");
  TCutG *b30ne = (TCutG*)fcut->Get("b30ne");
  TCutG *b31ne = (TCutG*)fcut->Get("b31ne");
  TCutG *b32ne = (TCutG*)fcut->Get("b32ne");
  TCutG *b29f = (TCutG*)fcut->Get("b29f");
  fcut->Close();
  
  
  TChain *chain = new TChain("FragPID","FragPID");
  chain->Add("./root/run029[5-9].root.FragPID");
  chain->Add("./root/run030[0-1].root.FragPID");

  
  TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";
  TCut trig_neut = "coinTrigger==3||coinTrigger==5||coinTrigger==7";
  TCut cut_f32ne = "abs(fragAoZ-3.34)<0.06&&abs(fragZ-3.95)<0.15";

  Double_t f32ne_x[9] = {3.32112,3.27802,3.26455,3.3292,3.35614,3.39116,3.40194,3.35075,3.32112};
  Double_t f32ne_y[9] = {4.04911,4.01403,3.90242,3.84503,3.83227,3.89286,4.00446,4.05548,4.04911};

  TCutG *f32ne = new TCutG("f32ne",9,f32ne_x,f32ne_y);
  f32ne->SetVarX("fragAoZ");
  f32ne->SetVarY("fragZ");

  
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  //c1->cd(1);
  //chain->Draw("beamZ:beamAoZ>>hist(400,2.8,3.3,400,7,13)","","colz");
  c1->cd(1);
  chain->Draw("fragZ:fragAoZ>>h1(400,2.0,3.5,400,3.5,5.0)","b34na"&&cut_reaction,"colz");
  h1->SetTitle("Fragment PID for{}^{34}Na beam w/o neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
  f32ne->Draw("SAME");
  c1->cd(2);
  chain->Draw("fragZ:fragAoZ>>h2(200,2.0,3.5,200,3.5,5.0)","b34na"&&cut_reaction&&trig_neut,"colz");
  h2->SetTitle("Fragment PID for{}^{34}Na beam w/ neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
  h2->GetZaxis()->SetRangeUser(1,100);

  /*  
      c1->cd(3);
      chain->Draw("fragZ:fragAoZ>>h3(400,2.0,3.5,400,3.5,5.0)","b34na"&&cut_reaction&&"f32ne","colz");
      h3->SetTitle("Fragment PID for{}^{34}Na beam w/o neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
      c1->cd(4);
      chain->Draw("fragZ:fragAoZ>>h4(200,2.0,3.5,200,3.5,5.0)","b34na"&&cut_reaction&&trig_neut&&"f32ne","colz");
      h4->SetTitle("Fragment PID for{}^{34}Na beam w/ neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
      h4->GetZaxis()->SetRangeUser(1,100);
  */
  
  /*
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->Divide(2,2);
    c1->cd(1);
    chain->Draw("fragZ:fragAoZ>>h1(400,2.0,3.5,400,3.5,5.0)","b34na","colz");
    c1->cd(2);
    chain->Draw("fragZ:fragAoZ>>h2(400,2.0,3.5,400,3.5,5.0)","b33na","colz");
    c1->cd(3);
    chain->Draw("fragZ:fragAoZ>>h3(400,2.0,3.5,400,3.5,5.0)","b32ne","colz");
    c1->cd(4);
    chain->Draw("fragZ:fragAoZ>>h4(400,2.0,3.5,400,3.5,5.0)","b31ne","colz");
  */
}
