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

  chain->Draw("f7PlaTL:sqrt(Q7)","","colz");

}
