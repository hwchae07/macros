void ped_nebula()
{
  TFile *file = new TFile("./root/run0303.root.NEBULA","r");
  TTree *tree = (TTree*)file->Get("NEBULA");

  std::ofstream fout("./dat/nebula_ped.dat");
  
  TH2 *h1;
  TH2 *h2;

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);

  c1->cd(1);
  tree->Draw("nebulaQURaw:nebulaID>>h1(144,0.5,144.5,400,0,400)","","colz");
  h1 = (TH2D*)gDirectory->Get("h1");
  c1->cd(2);
  tree->Draw("nebulaQDRaw:nebulaID>>h2(144,0.5,144.5,400,0,400)","","colz");
  h2 = (TH2D*)gDirectory->Get("h2");

  TH1 *hpro1;
  TH1 *hpro2;
  TF1 *fit1;
  TF1 *fit2;

  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);

  c2->Print("./fig/nebula_ped.pdf[");
  
  for(Int_t id=1;id<=144;id++)
    {
      c2->cd(1);
      hpro1 = (TH1D*)h1->ProjectionY(Form("ID%dQup",id),id,id);
      hpro1->SetTitle(Form("ID%d Qraw_{up};channel (ch)",id));
      hpro1->Fit("gaus","q0");
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
      hpro2 = (TH1D*)h2->ProjectionY(Form("ID%dQdown",id),id,id);
      hpro2->SetTitle(Form("ID%d Qraw_{down};channel (ch)",id));
      hpro2->Fit("gaus","q0");
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
      c2->Print("./fig/nebula_ped.pdf");
    }
  fout.close();
  c2->Print("./fig/nebula_ped.pdf]");
}
