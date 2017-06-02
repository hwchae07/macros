{
  TFile *plaTCut = new TFile("./cut/plaTCut.root","r");
  TCut *f3pla = f3plat;
  TCut *f5pla = f5plat;
  TCut *f7pla = f7plat;
  TCut *f13pla1 = f13plat1;
  TCut *f13pla2 = f13plat2;
  plaTCut->Close();
  
  TFile *file = new TFile("./root/run0275.root.BeamPla","r");
  TTree *tree = (TTree*)file->Get("BeamPla");
  tree->BuildIndex("RunNum","EventNum");
  //tree->AddFriend("BDC","./root/run0275.root.BDC");
  //tree->GetFriend("BDC")->BuildIndex("RunNum","EventNum");
  
  Double_t tof713_offset = 213.73055 + 352.345;   //nsec
  Double_t fl7_13 = (3957.1684 + 3957.2192) / 2.; //cm
  Double_t brho0 = 7.8412; //Tm

  tree->SetAlias("tofc713",Form("TOF7_13+%lf",tof713_offset));
  tree->SetAlias("betaF7",Form("%lf/tofc713/29.9792458",fl7_13));
  tree->SetAlias("betaF5","betaF7*1.005");
  tree->SetAlias("betaF13","betaF7*0.9");
  tree->SetAlias("dEFac_beam","TMath::Log(4866 * betaF13 * betaF13) - TMath::Log(1-betaF13*betaF13) - betaF13*betaF13");
  tree->SetAlias("beamZ","1+10*TMath::Sqrt(icbEloss/dEFac_beam)*betaF13");
  tree->SetAlias("beamMomU"," 931.494*betaF5/TMath::Sqrt(1-betaF5*betaF5)"); // mom/c = gamma*mass*beta [MeV/u/c]
  tree->SetAlias("beamAoZ",Form("(1+f5X/3300)*%lf * 299.792458 / beamMomU ",brho0));

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  tree->Draw("beamZ:beamAoZ>>(400,2.8,3.5,400,8,12)",f3pla[0]&&f5pla[0]&&f7pla[0]&&f13pla1[0]&&f13pla2[0],"colz");

}
