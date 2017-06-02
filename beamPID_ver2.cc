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
  chain->Add("./root/run027[5-9].root.BeamPID");
  //chain->Add("./root/run028[0-9].root.BeamPID");
  //chain->Add("./root/run029[0-4].root.BeamPID");

  //chain->Draw("beamZ:beamAoZ>>h1(1000,2.8,3.4,1000,6,12)","","colz");

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

}
