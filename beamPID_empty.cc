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
  
  TChain *chain = new TChain("BeamPID","BeamPID");
  chain->Add("./root/run029[5-9].root.BeamPID");
  chain->Add("./root/run030[0-1].root.BeamPID");

  TCut cut_reaction = "TMath::Sqrt(tgtX*tgtX + tgtY*tgtY) < 40";
  
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  chain->Draw("beamZ:beamAoZ>>h1(1000,2.8,3.4,1000,6,13)",cut_reaction,"colz");
  h1->SetTitle("Z:A/Z;A/Z;Z");

  
  TCanvas *c2 = new TCanvas("c2","c2",1200,400);
  c2->Divide(2,1);
  c2->cd(1);
  chain->Draw("beamAoZ>>h2(200,3.06,3.12)","b34na"&&cut_reaction);
  h2->Fit("gaus","rq","",3.08,3.095);
  h2->SetTitle("A/Z;A/Z");
  TF1 *fit1 = (TF1*)h2->GetFunction("gaus");
  cout<<"dA/A : "<<fit1->GetParameter(2)/fit1->GetParameter(1)*100<<endl;
  c2->cd(2);
  chain->Draw("beamZ>>h3(200,10.4,11.8)","b34na"&&cut_reaction);
  h3->SetTitle("Z;Z");
  h3->Fit("gaus","rq","",10.8,11.3);
  TF1 *fit2 = (TF1*)h3->GetFunction("gaus");
  cout<<"dZ/Z : "<<fit2->GetParameter(2)/fit2->GetParameter(1)*100<<endl;

  TCanvas *c3 = new TCanvas("c3","c3",1800,400);
  c3->Divide(3,1);
  c3->cd(1);
  chain->Draw("betaF7","b34na"&&cut_reaction);
  c3->cd(2);
  chain->Draw("beamMomU","b34na"&&cut_reaction);
  c3->cd(3);
  chain->Draw("beamKE","b34na"&&cut_reaction);
  
  //beam plastic slew effect check
  /*
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1);
    c1->cd(1);
    chain->Draw("f7PlaTL:sqrt(Q7)>>h1(200,14,30,200,658,677)","b34na","colz");
    TH1 *h1 = (TH1D*)gDirectory->Get("h1");
    h1->SetTitle("{}^{34}Na;#sqrt{Q};T (arbitrary)");
    c1->cd(2);
    chain->Draw("f7PlaTL:sqrt(Q7)>>h2(200,14,30,200,658,677)","b32ne","colz");
    TH1 *h2 = (TH1D*)gDirectory->Get("h2");
    h2->SetTitle("{}^{32}Ne;#sqrt{Q};T (arbitrary)");
    c1->Print("./fig/beamPlaSlewCheck.pdf");
  */
}
