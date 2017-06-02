{

  //beam cut 22Ne, 21F and fragment cut 22Ne//
  TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne, *b21f, *f22ne;
  f_cut->GetObject("b22ne",b22ne);
  f_cut->GetObject("b21f",b21f);
  //f_cut->GetObject("f22ne",f22ne);
  //beam cut 22Ne, 21F and fragment cut 22Ne//

  Int_t runNum = 133;
  TString filename = Form("run%04d.root",runNum);

  
  TChain *chain  = new TChain("BeamPla","BeamPla");
  chain->Add("./root/run013[3-8].root.BeamPla");
  TChain *chain_hod = new TChain("HOD","HOD");
  chain_hod->Add("./root/run013[3-8].root.HOD");
  //TChain *chain_fdc = new TChain("FDC","FDC");
  //chain_fdc->Add("./root/run013[3-8].root.FDC");

  
  chain->BuildIndex("RunNum","EventNum");
  chain->AddFriend(chain_hod);
  chain->GetFriend("HOD")->BuildIndex("RunNum","EventNum");
  //hain->AddFriend(chain_fdc);
  //chain->GetFriend("FDC")->BuildIndex("RunNum","EventNum");

  

  chain->SetAlias("fne","abs(hodQA-450)<30");
  for(Int_t i=1;i<=23;i++)
    {
      chain->SetAlias(Form("hod%02d_%02d",i,i+1),Form("hodNum==2&&((hodID[0]==%d&&hodID[1]==%d)||(hodID[0]==%d&&hodID[1]==%d))",i,i+1,i+1,i));
    }

  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  chain->Draw("hodTA>>h1(100,75,95)","b22ne&&hod11_12&&hodID==11");
  Double_t center1 = h1->GetBinCenter(h1->GetMaximumBin());
  //cout<<h1->GetMaximum()<<endl;
  //cout<<h1->GetRMS()<<endl;
  TF1 *fit1 = new TF1("fit1","gaus");
  fit1->SetParameter(0,h1->GetMaximum());
  fit1->SetParameter(1,center1);
  fit1->SetParameter(2,h1->GetRMS());
  h1->Fit("fit1","R","",center1-1,center1+1);

  c1->Update();
  
  
  c1->cd(2);
  chain->Draw("icbEloss:TOF3_13","b22ne","colz");
  /*
  chain->Draw("hodTA>>h2(100,75,95)","b22ne&&hod11_12&&hodID==12");
  Double_t center2 = h2->GetBinCenter(h2->GetMaximumBin());
  TF1 *fit2 = new TF1("fit2","gaus");
  fit2->SetParameter(0,h2->GetMaximum());
  fit2->SetParameter(1,center2);
  fit2->SetParameter(2,h2->GetRMS());
  h2->Fit("fit2","R","",center2-1,center2+1);
  */
  c1->Update();
  
  
  

}
