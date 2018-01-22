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

  
  TChain *chain = new TChain("HOD","HOD");
  //chain->Add("./root/run0275.root.HOD");
  chain->Add("./root/run027[5-9].root.HOD");
  chain->Add("./root/run028[0-9].root.HOD");
  chain->Add("./root/run029[0-4].root.HOD");


  
}
