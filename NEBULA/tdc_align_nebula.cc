{
  TFile *file = new TFile("./root/run0303.root.NEBULA","r");
  TTree *tree = (TTree*)file->Get("NEBULA");

  ofstream fout1("./dat/tdc_align_nebula_up.dat");
  ofstream fout2("./dat/tdc_align_nebula_down.dat");
  
  //TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  //c1->Divide(2,2);
  
  TH1 *h1;
  TH1 *h2;
  TH1 *h3;
  TH1 *h4;

  int id = 1;

  for(id = 1; id<=120;id++)
    {
      cout<<"ID"<<id<<"is ongoing..."<<endl;
      TString raw_tu = Form("nebulaTURaw>500&&nebulaTURaw<3500&&nebulaID==%d",id);
      TString raw_td = Form("nebulaTDRaw>500&&nebulaTDRaw<3500&&nebulaID==%d",id);

  
      //c1->cd(1);
      tree->Draw("nebulaTU>>h1(400,0,250)",raw_tu,"goff");
      h1 = (TH1D*)gDirectory->Get("h1");
      //c1->cd(3);
      tree->Draw("nebulaTU2>>h2(400,0,250)",raw_tu,"goff");
      h2 = (TH1D*)gDirectory->Get("h2");
  
      //c1->cd(2);
      tree->Draw("nebulaTD>>h3(400,0,250)",raw_td,"goff");
      h3 = (TH1D*)gDirectory->Get("h3");
      //c1->cd(4);
      tree->Draw("nebulaTD2>>h4(400,0,250)",raw_td,"goff");
      h4 = (TH1D*)gDirectory->Get("h4");

  
      fout1<< h2->GetMean() - h1->GetMean() << endl;
      fout2<< h4->GetMean() - h3->GetMean() << endl;
    }

  fout1.close();
  fout2.close();
  
}
