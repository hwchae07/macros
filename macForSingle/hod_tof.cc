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
  
  TFile *file = new TFile("./root/run0275.root.BeamPID","r");
  TTree *tree = (TTree*)file->Get("BeamPID");
  tree->BuildIndex("EventNum");
  
  tree->AddFriend("BDC","./root/run0275.root.BDC");
  tree->GetFriend("BDC")->BuildIndex("EventNum");
  tree->AddFriend("FDC","./root/run0275.root.FDC");
  tree->GetFriend("FDC")->BuildIndex("EventNum");
  tree->AddFriend("HOD","./root/run0275.root.HOD");
  tree->GetFriend("HOD")->BuildIndex("EventNum");
  
  tree->AddFriend("FragPID","./root/run0275.root.FragPID");
  tree->GetFriend("FragPID")->BuildIndex("EventNum");

  TCut cut_central = "abs(f5X)<10";
  
  //Target size//
  Double_t bdcWidth = 120.;
  Double_t bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
  Double_t bdc1Tgt = bdc2Tgt + 1000;
  tree->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
  tree->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
  TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";  
  //Target size//
  
  //FDC1 & FDC2//
  TCut cut_fdc1 = "sqrt(fdc1X**2 + fdc1Y**2) < 40";
  TCut cut_fdc2 = "((fdc2X-150)/100)**2 + (fdc2Y/200)**2 < 1";
  //FDC1 & FDC2//
  TCut cut_total = cut_central&&cut_reaction&&cut_fdc1&&cut_fdc2;

  TCut cut_f32ne = "abs(fragZ-3.85)<0.15&&abs(fragAoZ-3.38)<0.04";

  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  tree->Draw("(hodTU+hodTD)/2 - T13>>h1(100,-226,-224)",cut_total&&cut_f32ne&&"b32ne"&&"hodID==11");
  h1->GetXaxis()->SetNdivisions(505);
  h1->Fit("gaus","rq","",-225.05,-224.65);
  TF1 *fit1 = (TF1*)h1->GetFunction("gaus");
  cout<<endl;
  cout<<"Simulation result : 52.8"<<endl;
  cout<<"old result : "<<fit1->GetParameter(1)<<", offset -> "<<52.8-fit1->GetParameter(1)<<endl;
  
  
  c1->cd(2);
  tree->Draw("(hodTUCal+hodTDCal)/2 - T13>>h2",cut_total&&"b32ne"&&"hodID==11");
  
  tree->Draw("(hodTUCal+hodTDCal)/2 - T13>>h2(100,-252,-250)",cut_total&&cut_f32ne&&"b32ne"&&"hodID==11");
  h2->GetXaxis()->SetNdivisions(505);
  h2->Fit("gaus","rq","",-251.2,-250.8);
  TF1 *fit2 = (TF1*)h2->GetFunction("gaus");
  cout<<"new result : "<<fit2->GetParameter(1)<<", offset -> "<<52.8-fit2->GetParameter(1)<<endl;

  cout<<"difference : "<<fit1->GetParameter(1)-fit2->GetParameter(1)<<endl;
  
}
