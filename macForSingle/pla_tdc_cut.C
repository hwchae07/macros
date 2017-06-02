void pla_tdc_cut(Int_t detNum)
{
  gROOT->Reset();
  TFile *file = new TFile("./root/run0275.root.BeamPla","read");
  TTree *tree = (TTree*)file->Get("BeamPla");
  
  TString detNameList[5] = {"f3Pla","f5Pla","f7Pla","f13Pla1","f13Pla2"};
  TString detName;

  TCut plaTCutList1[5];
  TCut plaTCutList2[5];
  TCut plaTCut(Form("%splaTCut",detName.Data()),Form("%splaTCut",detName.Data()));

  detName = detNameList[detNum];

  plaTCutList1[0] = Form("abs((%sTLRaw-%sTRRaw)-20)<60",detName.Data(),detName.Data());
  plaTCutList1[1] = Form("abs((%sTLRaw-%sTRRaw)+200)<150",detName.Data(),detName.Data());
  plaTCutList1[2] = Form("abs((%sTLRaw-%sTRRaw)-200)<100",detName.Data(),detName.Data());
  plaTCutList1[3] = Form("abs((%sTLRaw-%sTRRaw))<100",detName.Data(),detName.Data());
  plaTCutList1[4] = Form("abs((%sTLRaw-%sTRRaw)+200)<100",detName.Data(),detName.Data());

  plaTCutList2[0] = Form("abs(%sTLRaw-45500)<1500&&abs(%sTRRaw-45500)<1500",detName.Data(),detName.Data());
  plaTCutList2[1] = Form("abs(%sTLRaw-42000)<1000&&abs(%sTRRaw-42000)<1000",detName.Data(),detName.Data());
  plaTCutList2[2] = Form("abs(%sTLRaw-42800)<600&&abs(%sTRRaw-42500)<800",detName.Data(),detName.Data());
  plaTCutList2[3] = Form("abs(%sTLRaw-19900)<100&&abs(%sTRRaw-19900)<100",detName.Data(),detName.Data());
  plaTCutList2[4] = Form("abs(%sTLRaw-19800)<200&&abs(%sTRRaw-19800)<200",detName.Data(),detName.Data());
  
  plaTCut = plaTCutList1[detNum]&&plaTCutList2[detNum];

  /*
  TFile *fcut = new TFile("./cut/plaTCut.root","update");
  plaTCut.Write();
  fcut->Close();
  */
  
  TCanvas *c1 = new TCanvas("c1","c1",1500,1000);
  c1->Divide(3,2);

  c1->cd(1);
  tree->Draw(Form("%sTLRaw:%sTRRaw>>h1",detName.Data(),detName.Data()),"","colz");
  TH2 *h1 = (TH2D*)gDirectory->Get("h1");
  h1->SetTitle(Form("%sTLRaw VS %sTRRaw;channel (ch);channel (ch)",detName.Data(),detName.Data()));
  h1->GetXaxis()->SetNdivisions(509);
  h1->SetTitleOffset(1.5,"Y");

  c1->cd(2);
  tree->Draw(Form("%sTLRaw-%sTRRaw>>h2",detName.Data(),detName.Data()),plaTCut);
  TH1 *h2 = (TH1D*)gDirectory->Get("h2");
  h2->SetTitle(Form("%sTLRaw - %sTRRaw;channel (ch)",detName.Data(),detName.Data()));

  c1->cd(4);
  tree->Draw(Form("%sTLRaw:%sTRRaw>>h3",detName.Data(),detName.Data()),plaTCut,"colz");
  TH2 *h3 = (TH2D*)gDirectory->Get("h3");
  h3->SetTitle(Form("%sTLRaw VS %sTRRaw;channel (ch);channel (ch)",detName.Data(),detName.Data()));
  h3->GetXaxis()->SetNdivisions(508);
  h3->SetTitleOffset(1.5,"Y");
  
  c1->cd(5);
  tree->Draw(Form("%sTLRaw>>h4",detName.Data()),plaTCut);
  TH1 *h4 = (TH1D*)gDirectory->Get("h4");
  h4->SetTitle(Form("%sTLRaw;channel (ch)",detName.Data()));
  h4->GetXaxis()->SetNdivisions(509);

  c1->cd(6);
  tree->Draw(Form("%sTRRaw>>h5",detName.Data()),plaTCut);
  TH1 *h5 = (TH1D*)gDirectory->Get("h5");
  h5->SetTitle(Form("%sTRRaw;channel (ch)",detName.Data()));
  h5->GetXaxis()->SetNdivisions(509);



  //c1->Print(Form("./fig/%sTDC_real.pdf",detName.Data()));
  
  
}


void gen_pla_tdc_cut()
{
  gROOT->Reset();
  TFile *file = new TFile("./root/run0275.root.BeamPla","read");
  TTree *tree = (TTree*)file->Get("BeamPla");
  
  TString detNameList[5] = {"f3Pla","f5Pla","f7Pla","f13Pla1","f13Pla2"};
  TString detName;

  TCut plaTCutList1[5];
  TCut plaTCutList2[5];
  TCut plaTCut[5];
  TCut f3plat[1];
  TCut f5plat[1];
  TCut f7plat[1];
  TCut f13plat1[1];
  TCut f13plat2[1];

  for(Int_t detNum=0;detNum<5;detNum++)
    {
  
      detName = detNameList[detNum];

      plaTCutList1[0] = Form("abs((%sTLRaw-%sTRRaw)-20)<60",detName.Data(),detName.Data());
      plaTCutList1[1] = Form("abs((%sTLRaw-%sTRRaw)+200)<150",detName.Data(),detName.Data());
      plaTCutList1[2] = Form("abs((%sTLRaw-%sTRRaw)-200)<100",detName.Data(),detName.Data());
      plaTCutList1[3] = Form("abs((%sTLRaw-%sTRRaw))<100",detName.Data(),detName.Data());
      plaTCutList1[4] = Form("abs((%sTLRaw-%sTRRaw)+200)<100",detName.Data(),detName.Data());

      plaTCutList2[0] = Form("abs(%sTLRaw-45500)<1500&&abs(%sTRRaw-45500)<1500",detName.Data(),detName.Data());
      plaTCutList2[1] = Form("abs(%sTLRaw-42000)<1000&&abs(%sTRRaw-42000)<1000",detName.Data(),detName.Data());
      plaTCutList2[2] = Form("abs(%sTLRaw-42800)<600&&abs(%sTRRaw-42500)<800",detName.Data(),detName.Data());
      plaTCutList2[3] = Form("abs(%sTLRaw-19900)<100&&abs(%sTRRaw-19900)<100",detName.Data(),detName.Data());
      plaTCutList2[4] = Form("abs(%sTLRaw-19800)<200&&abs(%sTRRaw-19800)<200",detName.Data(),detName.Data());
  
      plaTCut[detNum] = plaTCutList1[detNum]&&plaTCutList2[detNum];
      //plaTCut[detNum].SetTitle(Form("%sTDC_cut",detName.Data()));
      plaTCut[detNum].SetName(Form("%sTDC_cut",detName.Data()));
    }

  
  f3plat[0] = plaTCut[0];
  f5plat[0] = plaTCut[1];
  f7plat[0] = plaTCut[2];
  f13plat1[0] = plaTCut[3];
  f13plat2[0] = plaTCut[4];

  f3plat[0].SetName("f3plat");
  f5plat[0].SetName("f5plat");
  f7plat[0].SetName("f7plat");
  f13plat1[0].SetName("f13plat1");
  f13plat2[0].SetName("f13plat2");
  
  TFile *fcut = new TFile("./cut/plaTCut.root","recreate");
  fcut->WriteTObject(f3plat);
  fcut->WriteTObject(f5plat);
  fcut->WriteTObject(f7plat);
  fcut->WriteTObject(f13plat1);
  fcut->WriteTObject(f13plat2);
  fcut->Close();
  
}

