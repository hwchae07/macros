{
  TChain *chain1 = new TChain("HOD","HOD");
  TChain *chain2 = new TChain("BeamPla","BeamPla");
  
  chain1->Add("./root/run013[3-4].root.HOD");
  chain2->Add("./root/run013[3-4].root.BeamPla");

  //  TFile *file1 = new TFile("root/run0133.root.HOD");
  //  TTree *chain1 = (TTree*)file1->Get("HOD");
  //  TFile *file2 = new TFile("root/run0133.root.BeamPla");
  //  TTree *chain2 = (TTree*)file2->Get("BeamPla");
  
  chain1->BuildIndex("RunNum","EventNum");
  //  chain1->GetFriend("BeamPla")->BuildIndex("RunNum","EventNum");
  chain2->BuildIndex("RunNum","EventNum");
  chain1->AddFriend(chain2);

  Int_t runNum, eventNum;
  Int_t runNum2, eventNum2;
  Int_t hodNum;
  chain1->SetBranchAddress("RunNum",&runNum);
  chain1->SetBranchAddress("EventNum",&eventNum);
  chain1->SetBranchAddress("hodNum",&hodNum);
  chain1->SetBranchAddress("hodID",&hodID);
  chain1->SetBranchAddress("hodTA",&hodTA);

  
  ifstream file("./dat/hod_offset.dat");

  Int_t id[25];
  Double_t offset[25];
  for(Int_t i=0;i<=24;i++)
    {
      file>>id[i]>>offset[i];
      //cout<<id[i]<<"\t"<<offset[i]<<endl;
    }

  //  for(Int_t i=0;i<chain1->GetEntries();i++)
  for(Int_t i=0;i<1000;i++)
    {
      if (!(i % 100)){
	cout << "\rEvent : " << i << " / Progress: " << (Double_t)i/chain1->GetEntries()*100. << "%";
	cout.flush();}
      chain1->GetEntry(i);
      //chain2->GetEntry(i);

      //////////////////////////////////////////////////
      // Check all friends
      Bool_t goodFlag = true;
      TIter nextf(chain1->GetListOfFriends());
      TFriendElement* fe = 0;
      while ((fe = (TFriendElement*) nextf())) {
	TTree *t = fe->GetTree();
	if (t && t->GetEntryNumberWithIndex(runNum,eventNum) < 0)
	  goodFlag = false;}
      
      if (!goodFlag || chain1->GetEntryWithIndex(runNum,eventNum) <= 0) 
	continue;
      //      chain2->GetEntryWithIndex("RunNum","EventNum"));
      
      if (eventNum != eventNum2){ 
	cout << runNum << " ";
	cout << runNum2 << endl;
	cout << eventNum << " ";
	cout << eventNum2 << endl;
	cout << hodNum << endl;}
      //chain1->Scan("EventNum");
      //chain2->Scan("EventNum");
      
    }
  
  file.close();
}
