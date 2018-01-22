{

  //beam cut 22Ne, 21F//
  TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne, *b21f;
  f_cut->GetObject("b22ne",b22ne);
  f_cut->GetObject("b21f",b21f);
  //beam cut 22Ne, 21F//
  
  Int_t runNum = 133;
  TString filename = Form("run%04d.root",runNum);

  
  TChain *chain  = new TChain("BeamPla","BeamPla");
  chain->Add("./root/run013[3-4].root.BeamPla");
  chain->BuildIndex("RunNum","EventNum");
  TChain *chain_hod  = new TChain("HOD","HOD");
  chain_hod->Add("./root/run013[3-4].root.HOD");
  chain_hod->BuildIndex("RunNum","EventNum");
  

  chain->AddFriend(chain_hod);
  //  chain->GetFriend("HOD")->BuildIndex("RunNum","EventNum");
  //chain->AddFriend("HOD","./root/run0133.root.HOD");
  //chain->GetFriend("HOD")->BuildIndex("RunNum","EventNum");

  //chain->Draw("icbEloss:TOF3_13>>(200,-440,-380,200,0,50)","b22ne","colz");
  //  chain->Draw("hodTA:hodID>>(24,0.5,24.5,200,50,100)","b22ne&&RunNum==133","colz");

  
  

}
