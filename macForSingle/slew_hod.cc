{
  TFile *file = new TFile("./root/run0275.root.BeamPla","r");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("EventNum");
  tree->AddFriend("HOD","./root/run0275.root.HOD","r");
  tree->GetFriend("HOD")->BuildIndex("EventNum");
}
