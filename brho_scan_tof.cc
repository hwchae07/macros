{

  //load the hod offset//
  ifstream file("./dat/hod_offset_brho_scan.dat");
  //ifstream file("./dat/hod_offset.dat");

  Int_t dummy[25];
  Double_t offset[25];
  for(Int_t i=0;i<=24;i++)
    file>>dummy[i]>>offset[i];

  TString tof_hod = "";

  tof_hod += "hodTA-T13";
  /*
  tof_hod += "+";
  tof_hod += offset[0];

  for(Int_t id=1;id<=24;id++)
    {
      tof_hod += "+";
      tof_hod += offset[id];
      tof_hod += Form("*(hodID==%d)",id);
    }
  */
  //load the hod offset//


  
  TFile *file3 = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne = (TCutG*)file3->Get("b22ne");
  TCut cut_central = "abs(f5X)<2";
  
  TH2 *hist1 = new TH2D("hist1","hist1",24,0.5,24.5,200,30,80);
  TH2 *hsum1 = new TH2D("hsum1","hsum1",24,0.5,24.5,200,30,80);

  TH2 *hist2 = new TH2D("hist2","hist2",24,0.5,24.5,1000,-250,-200);
  TH2 *hsum2 = new TH2D("hsum2","hsum2",24,0.5,24.5,1000,-250,-200);
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  for(Int_t runNum = 133;runNum<=138;runNum++)
    {
    
      TFile *file1 = new TFile(Form("./root/run%04d.root.HOD",runNum));
      TTree *tree1 = (TTree*)file1->Get("HOD");
      TFile *file2 = new TFile(Form("./root/run%04d.root.BeamPla",runNum));
      TTree *tree2 = (TTree*)file2->Get("BeamPla");
  
      tree1->BuildIndex("RunNum","EventNum");
      tree1->AddFriend(tree2);
      tree1->GetFriend("BeamPla")->BuildIndex("RunNum","EventNum");

      tree1->SetAlias("tof_hod",tof_hod);

      tree1->Draw("tof_hod:hodID>>hist1(24,0.5,24.5,200,30,80)","b22ne","colz");
      hist1=(TH2D*)gDirectory->Get("hist1");
      hsum1->Add(hist1);
      hsum1->Draw("colz");
      
      /*    
	    tree1->Draw("hodTA-T13:hodID>>hist2(24,0.5,24.5,1000,-250,-200)","b22ne"&&cut_central,"colz");
	    hist2=(TH2D*)gDirectory->Get("hist2");
	    hsum2->Add(hist2);
	    hsum2->Draw("colz");
      */
      c1->Update();
    }
  
  
  /*
  //save the histogram//
  TFile *file2 = new TFile("./hist/brho_scan_hodt.root","recreate");
  hsum2->Write();
  file2->Close();
  //save the histogram//
  */

  
  /*
  //fitting process//
  TFile *file4 = new TFile("./hist/brho_scan_hodt.root","r");
  TH2 *h1;
  TH1 *hpy;
  TF1 *fit1;
  file4->GetObject("hsum2",h1);

  TCanvas *c2 = new TCanvas("c2","c2",1800,1200);
  c2->Divide(6,4);

  for(Int_t id=1;id<=24;id++)
  {
  c2->cd(id);
  h1->ProjectionY(Form("ID%d",id),id,id);
  gDirectory->GetObject(Form("ID%d",id),hpy);
  hpy->Rebin(4);
  hpy->SetTitle(Form("ID%d;TOF",id));
  hpy->Fit("gaus","q","");
  hpy->GetXaxis()->SetRangeUser(-230,-220);
  hpy->Draw();
  fit1 = hpy->GetFunction("gaus");
  cout<<"ID"<<id<<"\t"<<fit1->GetParameter(1)<<endl;
  }
  //fitting process//
  */
  file.close();
}
