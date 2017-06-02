{

  //beam cut 22Ne, 21F and fragment cut 22Ne//
  TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne, *b21f, *f22ne;
  f_cut->GetObject("b22ne",b22ne);
  f_cut->GetObject("b21f",b21f);
  f_cut->GetObject("f22ne",f22ne);
  //beam cut 22Ne, 21F and fragment cut 22Ne//



  
  Int_t runNum = 133;
  TString filename = Form("run%04d.root",runNum);

  TFile *file = new TFile("./root/run0133.root.BeamPla","R");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("RunNum","EventNum");
  tree->AddFriend("FDC","./root/run0133.root.FDC");
  tree->GetFriend("FDC")->BuildIndex("RunNum","EventNum");
  tree->AddFriend("HOD","./root/run0133.root.HOD");
  tree->GetFriend("HOD")->BuildIndex("RunNum","EventNum");

  tree->Draw("hodQA:fdc2X>>h1(200,-1200,-600,200,200,500)","b22ne&&f22ne","colz");

}
