void ped()
{
  TFile *file = new TFile("./root/run0302.root.HOD","r");
  TTree *tree = (TTree*)file->Get("HOD");
  std::ofstream fout("./dat/hod_ped.dat");

  TH2 *h1;
  TH2 *h2;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  tree->Draw("hodQURaw:hodID>>h1(24,0.5,24.5,400,0,400)","","colz");
  h1 = (TH2D*)gDirectory->Get("h1");
  c1->cd(2);
  tree->Draw("hodQDRaw:hodID>>h2(24,0.5,24.5,400,0,400)","","colz");
  h2 = (TH2D*)gDirectory->Get("h2");

  TH1 *hpro1;
  TH1 *hpro2;
  TF1 *fit1;
  TF1 *fit2;

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  for(int i=1;i<=24;i++)
    {
      c2->cd(1);
      hpro1 = (TH1D*)h1->ProjectionY(Form("ID%dQU",i),i,i);
      hpro1->Fit("gaus","q");
      fit1 = (TF1*)hpro1->GetFunction("gaus");
      Double_t sup1 = fit1->GetParameter(1) + 10*fit1->GetParameter(2);
      Double_t inf1 = fit1->GetParameter(1) - 10*fit1->GetParameter(2);
      Double_t Rrange1 = fit1->GetParameter(1) + 3*fit1->GetParameter(2);
      Double_t Lrange1 = fit1->GetParameter(1) - 3*fit1->GetParameter(2);
      hpro1->GetXaxis()->SetRangeUser(inf1,sup1);
      hpro1->Fit("gaus","rq","",Lrange1,Rrange1);
      fit1 = (TF1*)hpro1->GetFunction("gaus");
      hpro1->Draw();
      c2->cd(2);
      hpro2 = (TH1D*)h2->ProjectionY(Form("ID%dQD",i),i,i);
      hpro2->Fit("gaus","q");
      fit2 = (TF1*)hpro2->GetFunction("gaus");
      Double_t sup2 = fit2->GetParameter(1) + 10*fit2->GetParameter(2);
      Double_t inf2 = fit2->GetParameter(1) - 10*fit2->GetParameter(2);
      Double_t Rrange2 = fit2->GetParameter(1) + 3*fit2->GetParameter(2);
      Double_t Lrange2 = fit2->GetParameter(1) - 3*fit2->GetParameter(2);
      hpro2->GetXaxis()->SetRangeUser(inf2,sup2);
      hpro2->Fit("gaus","rq","",Lrange2,Rrange2);
      fit2 = (TF1*)hpro2->GetFunction("gaus");
      hpro2->Draw();
      c2->Update();

      fout<<setw(6)<<fit1->GetParameter(1)<<" "<<fit2->GetParameter(1)<<endl;
      
      getchar();
    }

  fout.close();
}

void check_hod()
{
  TFile *file = new TFile("./root/run0302.root.HOD","r");
  TTree *tree = (TTree*)file->Get("HOD");

  TH2 *h1;
  TH2 *h2;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->cd(1);
  tree->Draw("hodQURaw:hodID>>h1(24,0.5,24.5,400,0,400)","","colz");
  h1 = (TH2D*)gDirectory->Get("h1");
  h1->SetTitle("HODF24 QDC_{U,Raw} VS ID;ID;QDC (channel)");
  c1->cd(2);
  tree->Draw("hodQU:hodID>>h2(24,0.5,24.5,400,-100,300)","","colz");
  h2 = (TH2D*)gDirectory->Get("h2");
  h2->SetTitle("HODF24 QDC_{U,Ped} VS ID;ID;QDC (channel)");
}
