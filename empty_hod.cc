{
  Int_t runNum = 295;
  TString filename = Form("run%04d.root",runNum);
  
  TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("EventNum");
  tree->AddFriend("BDC",Form("./root/%s.BDC",filename.Data()));
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
  tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  
  //Beam PID//
  tree->SetAlias("tofc","TOF3_13-3*dT5");
  tree->SetAlias("elossc","icbEloss-0.1*tofc-40");
  TCut cut_na = "TMath::Abs(elossc-22.5)<2.5";
  TCut cut_34na = cut_na&&"TMath::Abs(tofc+388)<3";
  TCut cut_ne = "abs(elossc-17.8)<2.2";
  TCut cut_32ne = cut_ne&&"abs(tofc+379)<3";
  //Beam PID//

  TCut cut_central = "abs(f5X)<10";

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  
  tree->Draw("hodTA - T13>>h1(100,-228,-225)",cut_32ne&&cut_central&&"hodID==8");
  TH1 *h1;
  gDirectory->GetObject("h1",h1);
  h1->Fit("gaus","qr","",-226.7,-226.3);


  //HOD time offset//
  /*ifstream file1("./dat/hod_offset_empty.dat");

  Int_t dummy[25];
  Double_t offset[25];
  for(Int_t i=0;i<=24;i++)
    file1>>dummy[i]>>offset[i];
  */
  TString tof_hod = "";

  tof_hod += "hodTA-T13";
  /*tof_hod += "+";
  tof_hod += offset[0];

  for(Int_t id=1;id<=24;id++)
    {
      tof_hod += "+";
      tof_hod += offset[id];
      tof_hod += Form("*(hodID==%d)",id);
    }
  */  
  tree->SetAlias("tof_hod",tof_hod);
  //HOD time offset//
  
  tree->Draw("tof_hod:fdc2X>>(200,0,500,200,45,55)",cut_32ne&&"fdc2X>0&&abs(tof_hod-50)<10","colz");

  //  tree->Draw("fdc2X:hodID>>(200,0,30,200,0,500)",cut_32ne,"colz");
  
      
}
