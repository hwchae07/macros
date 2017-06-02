void BF(Int_t runNum,
	const char *det1,
	const char *det2 = 0,
	const char *det3=0,
	const char *det4=0,
	const char *det5=0,
	const char *det6=0,
	const char *det7=0,
	const char *det8=0,
	const char *det9=0,
	const char *det10=0)
{
  TFile *file = new TFile(Form("./root/run%04d.root.%s",runNum,det1),"read");
  TTree *tree = (TTree*)file->Get(Form("%s",det1));
  tree->BuildIndex("RunNum","EventNum");

  vector<const char*> dets;
  
  if(det2)  dets.push_back(det2);
  if(det3)  dets.push_back(det3);
  if(det4)  dets.push_back(det4);
  if(det5)  dets.push_back(det5);
  if(det6)  dets.push_back(det6);
  if(det7)  dets.push_back(det7);
  if(det8)  dets.push_back(det8);
  if(det9)  dets.push_back(det9);
  if(det10) dets.push_back(det10);

  for(Int_t i=0;i<dets.size();i++)
    {
      tree->AddFriend(dets[i],Form("./root/run%04d.root.%s",runNum,dets[i]));
      tree->GetFriend(dets[i])->BuildIndex("RunNum","EventNum");
    }
}
