//after position calibration we can do timing correction using position information
//using Q after pedestal calibration

//function?
//void save_hist() : make 2d histogram of tof vs Q for each ID

TCutG* makeCUT(TF1 *func, double width, int pointNum, double *cutX, double *cutY)
{
  cutX = new double[pointNum];
  cutY = new double[pointNum];

  for(int i=0; i<pointNum/2 ; i++)
    {
      cutX[i]    = 10 + i*(8000./pointNum);
      cutY[i]    = func->Eval(cutX[i]) + width;

      cutX[i+pointNum/2] = 10 + 4000 - (i+1)*(8000./pointNum);
      double temp = 10 + 4000 - (i+1)*(8000./pointNum);
      cutY[i+pointNum/2] = func->Eval(temp) - width;
    }

  TCutG *cut = new TCutG("cut",pointNum,cutX,cutY);
  return cut;
}

void slew_res_check()
{
  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  TH1 *hu1;
  TH1 *hu2;
  TH1 *hu3;
  TH1 *hu4;
  TH1 *hd1;
  TH1 *hd2;
  TH1 *hd3;
  TH1 *hd4;

  TCanvas *c1 = new TCanvas("c1","c1",1800,900);
  c1->Divide(4,2);

  c1->cd(1);
  chain->Draw("nebulaTOFUGPS>>hu1(200,-5,5)","nebulaID>=1&&nebulaID<=30&&nebulaQA>6");
  hu1 = (TH1D*)gDirectory->Get("hu1");
  hu1->Fit("gaus","rq");
  fitu1 = (TF1*)hu1->GetFunction("gaus");
  hu1->Fit("gaus","rq","",fitu1->GetParameter(1)-2*fitu1->GetParameter(2),fitu1->GetParameter(1)+2*fitu1->GetParameter(2));
  hu1->SetTitle("NEBULA layer 1 (UP);TOF_{#gamma} (nsec);");

  c1->cd(2);
  chain->Draw("nebulaTOFUGPS>>hu2(200,-5,5)","nebulaID>=30&&nebulaID<=60&&nebulaQA>6");
  hu2 = (TH1D*)gDirectory->Get("hu2");
  hu2->Fit("gaus","rq");
  fitu2 = (TF1*)hu2->GetFunction("gaus");
  hu2->Fit("gaus","rq","",fitu2->GetParameter(1)-2*fitu2->GetParameter(2),fitu2->GetParameter(1)+2*fitu2->GetParameter(2));
  hu2->SetTitle("NEBULA layer 2 (UP);TOF_{#gamma} (nsec);");

  c1->cd(3);
  chain->Draw("nebulaTOFUGPS>>hu3(200,-5,5)","nebulaID>=60&&nebulaID<=90&&nebulaQA>6");
  hu3 = (TH1D*)gDirectory->Get("hu3");
  hu3->Fit("gaus","rq");
  fitu3 = (TF1*)hu3->GetFunction("gaus");
  hu3->Fit("gaus","rq","",fitu3->GetParameter(1)-2*fitu3->GetParameter(2),fitu3->GetParameter(1)+2*fitu3->GetParameter(2));
  hu3->SetTitle("NEBULA layer 3 (UP);TOF_{#gamma} (nsec);");

  c1->cd(4);
  chain->Draw("nebulaTOFUGPS>>hu4(200,-5,5)","nebulaID>=90&&nebulaID<=120&&nebulaQA>6");
  hu4 = (TH1D*)gDirectory->Get("hu4");
  hu4->Fit("gaus","rq");
  fitu4 = (TF1*)hu4->GetFunction("gaus");
  hu4->Fit("gaus","rq","",fitu4->GetParameter(1)-2*fitu4->GetParameter(2),fitu4->GetParameter(1)+2*fitu4->GetParameter(2));
  hu4->SetTitle("NEBULA layer 4 (UP);TOF_{#gamma} (nsec);");


  c1->cd(5);
  chain->Draw("nebulaTOFDGPS>>hd1(200,-5,5)","nebulaID>=1&&nebulaID<=30&&nebulaQA>6");
  hd1 = (TH1D*)gDirectory->Get("hd1");
  hd1->Fit("gaus","rq");
  fitd1 = (TF1*)hd1->GetFunction("gaus");
  hd1->Fit("gaus","rq","",fitd1->GetParameter(1)-2*fitd1->GetParameter(2),fitd1->GetParameter(1)+2*fitd1->GetParameter(2));
  hd1->SetTitle("NEBULA layer 1 (DOWN);TOF_{#gamma} (nsec);");

  c1->cd(6);
  chain->Draw("nebulaTOFDGPS>>hd2(200,-5,5)","nebulaID>=30&&nebulaID<=60&&nebulaQA>6");
  hd2 = (TH1D*)gDirectory->Get("hd2");
  hd2->Fit("gaus","rq");
  fitd2 = (TF1*)hd2->GetFunction("gaus");
  hd2->Fit("gaus","rq","",fitd2->GetParameter(1)-2*fitd2->GetParameter(2),fitd2->GetParameter(1)+2*fitd2->GetParameter(2));
  hd2->SetTitle("NEBULA layer 2 (DOWN);TOF_{#gamma} (nsec);");

  c1->cd(7);
  chain->Draw("nebulaTOFDGPS>>hd3(200,-5,5)","nebulaID>=60&&nebulaID<=90&&nebulaQA>6");
  hd3 = (TH1D*)gDirectory->Get("hd3");
  hd3->Fit("gaus","rq");
  fitd3 = (TF1*)hd3->GetFunction("gaus");
  hd3->Fit("gaus","rq","",fitd3->GetParameter(1)-2*fitd3->GetParameter(2),fitd3->GetParameter(1)+2*fitd3->GetParameter(2));
  hd3->SetTitle("NEBULA layer 3 (DOWN);TOF_{#gamma} (nsec);");

  c1->cd(8);
  chain->Draw("nebulaTOFDGPS>>hd4(200,-5,5)","nebulaID>=90&&nebulaID<=120&&nebulaQA>6");
  hd4 = (TH1D*)gDirectory->Get("hd4");
  hd4->Fit("gaus","rq");
  fitd4 = (TF1*)hd4->GetFunction("gaus");
  hd4->Fit("gaus","rq","",fitd4->GetParameter(1)-2*fitd4->GetParameter(2),fitd4->GetParameter(1)+2*fitd4->GetParameter(2));
  hd4->SetTitle("NEBULA layer 4 (DOWN);TOF_{#gamma} (nsec);");

  c1->Print("./fig/NEBULA/slew_resoultion_layer.pdf");

}

void slew_check()
{
  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  ifstream fin1("./dat/slew_2nd_nebula_up_pol1_layer0.dat");
  ifstream fin2("./dat/slew_2nd_nebula_down_pol1_layer0.dat");


  double par0u[120]={0,};
  double par1u[120]={0,};
  double par0d[120]={0,};
  double par1d[120]={0,};

  for(int i = 0;i<120;i++)
    {
      fin1>>par0u[i]>>par1u[i];
      fin2>>par0d[i]>>par1d[i];
    }
  fin1.close();
  fin2.close();

  TH2 *hu;
  TH2 *hd;
  TH2 *huu;
  TH2 *hdd;
  TF1 *fitu = new TF1("fitu","[0] + [1]*TMath::Power(x,-0.5)",10,4000);
  TF1 *fitd = new TF1("fitd","[0] + [1]*TMath::Power(x,-0.5)",10,4000);

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2);


  int id = 15;

  c1->Print("./fig/NEBULA/slew_2nd_nebula_pol1_check.pdf[");

  for(id=1;id<=120;id++)
    {
      c1->cd(1);
      chain->Draw(Form("nebulaTOFUGP:nebulaQPed>>hu%d(1000,0,4000,200,-20,40)",id),Form("nebulaID==%d",id),"colz");
      hu = (TH2D*)gDirectory->Get(Form("hu%d",id));
      hu->SetTitle(Form("ID%d (Up)",id));
      hu->GetZaxis()->SetRangeUser(0.01,10);
      fitu->SetParameter(0,par0u[id-1]);
      fitu->SetParameter(1,par1u[id-1]);
      fitu->SetLineColor(4);
      fitu->SetLineStyle(2);
      fitu->SetLineWidth(1);
      fitu->Draw("SAME");

      c1->cd(2);
      chain->Draw(Form("nebulaTOFDGP:nebulaQPed>>hd%d(1000,0,4000,200,-20,40)",id),Form("nebulaID==%d",id),"colz");
      hd = (TH2D*)gDirectory->Get(Form("hd%d",id));
      hd->SetTitle(Form("ID%d (Down)",id));
      hd->GetZaxis()->SetRangeUser(0.01,10);
      fitd->SetParameter(0,par0d[id-1]);
      fitd->SetParameter(1,par1d[id-1]);
      fitd->SetLineColor(4);
      fitd->SetLineStyle(2);
      fitd->SetLineWidth(1);
      fitd->Draw("SAME");

      c1->cd(3);
      chain->Draw(Form("nebulaTOFUGPS:nebulaQPed>>huu%d(1000,0,4000,200,-20,40)",id),Form("nebulaID==%d",id),"colz");
      huu = (TH2D*)gDirectory->Get(Form("huu%d",id));
      huu->SetTitle(Form("ID%d (Up)",id));
      huu->GetZaxis()->SetRangeUser(0.01,10);

      c1->cd(4);

      chain->Draw(Form("nebulaTOFDGPS:nebulaQPed>>hdd%d(1000,0,4000,200,-20,40)",id),Form("nebulaID==%d",id),"colz");
      hdd = (TH2D*)gDirectory->Get(Form("hdd%d",id));
      hdd->SetTitle(Form("ID%d (Down)",id));
      hdd->GetZaxis()->SetRangeUser(0.01,10);

      c1->Update();
      c1->Print("./fig/NEBULA/slew_2nd_nebula_pol1_check.pdf");
    }
  c1->Print("./fig/NEBULA/slew_2nd_nebula_pol1_check.pdf]");
}

void slew_comp()
{
  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2);

  TH2* h1;
  TH2* h2;
  TH2* h3;
  TH2* h4;

  int id = 26;
  c1->cd(1);
  chain->Draw("nebulaTOFUG:nebulaQPed>>h1(1000,0,3000,100,-20,80)",Form("nebulaTOFU>-9000&&nebulaID==%d",id),"colz");
  h1 = (TH2D*)gDirectory->Get("h1");
  h1->GetZaxis()->SetRangeUser(0.1,10);
  h1->SetTitle(Form("TOF_{#gamma} VS Q_{ped} for ID %d;Q_{ped};TOF_{#gamma}",id));
  c1->cd(2);
  chain->Draw("nebulaTOFUG:nebulaQPed>>h2(1000,0,3000,100,-20,80)",Form("nebulaTOFU>-9000&&nebulaID==%d&&abs(nebulaY)<60",id),"colz");
  h2 = (TH2D*)gDirectory->Get("h2");
  h2->GetZaxis()->SetRangeUser(0.1,10);
  h2->SetTitle(Form("TOF_{#gamma} VS Q_{ped} for ID %d, |Y| < 60 mm;Q_{ped};TOF_{#gamma}",id));
  c1->cd(3);
  chain->Draw("nebulaTOFUG:nebulaQPed>>h3(1000,0,3000,100,-20,80)","nebulaTOFU>-9000&&nebulaID<31&&abs(nebulaY)<60","colz");
  h3 = (TH2D*)gDirectory->Get("h3");
  h3->GetZaxis()->SetRangeUser(0.1,10);
  h3->SetTitle("TOF_{#gamma} VS Q_{ped} for Layer 1, |Y| < 60 mm;Q_{ped};TOF_{#gamma}");
  c1->cd(4);
  chain->Draw("nebulaTOFUGP:nebulaQPed>>h4(1000,0,3000,100,-20,80)",Form("nebulaTOFU>-9000&&nebulaID==%d",id),"colz");
  h4 = (TH2D*)gDirectory->Get("h4");
  h4->GetZaxis()->SetRangeUser(0.1,10);
  h4->SetTitle(Form("TOF_{#gamma,y_{cor}} VS Q_{ped} for ID %d;Q_{ped};TOF_{#gamma,y_{cor}}",id));

}


void slew_cor(int layer=0)
{
  TFile *file = new TFile("./hist/NEBULA/slew_2nd_nebula.root","r");

  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  ofstream fout1(Form("./dat/slew_2nd_nebula_up_pol1_layer%d.dat",layer));
  ofstream fout2(Form("./dat/slew_2nd_nebula_down_pol1_layer%d.dat",layer));

  TH2 *h1;
  TH2 *h2;
  TH2 *h3;

  TH1 *h2p;
  TH1 *h3p;

  TH2 *hu;
  TH2 *hd;

  TH2 *huu;
  TH2 *hdd;

  TH2 *huuu;
  TH2 *hddd;

  TProfile *hpu;
  TProfile *hpd;

  TProfile *hpuu;
  TProfile *hpdd;

  TProfile *hpuuu;
  TProfile *hpddd;

  TSpectrum *ts1 = new TSpectrum();
  TSpectrum *ts2 = new TSpectrum();

  TF1 *f1;
  TF1 *f2;

  int id = 15;

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  TCanvas *c2 = new TCanvas("c2","c2",1200,800);
  c2->Divide(3,2);

  c1->Print(Form("./fig/NEBULA/slew_2nd_nebula_pol1_layer%d_1.pdf[",layer));
  c2->Print(Form("./fig/NEBULA/slew_2nd_nebula_pol1_layer%d_2.pdf[",layer));

  int startID = 1;
  int endID = 120;

  if(layer==1)
    {
      startID=1;
      endID=30;
    }
  else if(layer==2)
    {
      startID=31;
      endID=60;
    }
  else if(layer==3)
    {
      startID=61;
      endID=90;
    }
  else if(layer==4)
    {
      startID=91;
      endID=120;
    }

  for(id = startID; id <= endID ; id+=1)
    {
      cout<<"ID"<<id<< " is ongoing..."<<"\r";
      cout.flush();

      h1 = (TH2D*)file->Get(Form("h1id%d",id));
      h2 = (TH2D*)file->Get(Form("h2id%d",id));
      h3 = (TH2D*)file->Get(Form("h3id%d",id));

      c1->cd(1);
      h2p = (TH1D*)h2->ProjectionY(Form("h2p%d",id),0,4000);
      h2p->Rebin(10);
      ts1->Search(h2p);
      double center1 = ts1->GetPositionX()[0];
      h2p->Fit("gaus","rq0","goff",center1-10,center1+10);
      f1 = (TF1*)h2p->GetFunction("gaus");
      double urange = f1->GetParameter(1)-3.5*f1->GetParameter(2);
      h2p->Draw();

      c1->cd(2);
      h3p = (TH1D*)h3->ProjectionY(Form("h3p%d",id),0,4000);
      h3p->Rebin(10);
      ts2->Search(h3p);
      double center2 = ts2->GetPositionX()[0];
      h3p->Fit("gaus","rq","goff",center2-10,center2+10);
      f2 = (TF1*)h3p->GetFunction("gaus");
      double drange = f2->GetParameter(1)-3.5*f2->GetParameter(2);
      h3p->Draw();



      c2->cd(1);
      chain->Draw(Form("nebulaTOFUGP:nebulaQPed>>hu%d(100,0,4000,500,-20,80)",id),Form("nebulaID==%d&&abs(nebulaTOFUGP-%lf+8)<7",id,urange),"colz");
      hu = (TH2D*)gDirectory->Get(Form("hu%d",id));
      hpu = (TProfile*)hu->ProfileX(Form("hpu%d",id));
      TF1 *fitu = new TF1("fitu","[0]+[1]*TMath::Power(x,-0.5)",0,4000);
      fitu->SetParameter(1,40);
      fitu->SetParLimits(1,20,60);
      hpu->Fit("fitu","rq","same",10,500);

      hu->GetYaxis()->SetRangeUser(urange-20,urange+5);
      hu->GetZaxis()->SetRangeUser(0.01,100);


      //TCutG* makeCUT(TF1 *func, double width, int pointNum, double *cutX, double *cutY)

      double *cutXU;
      double *cutYU;
      TCutG* cutU = makeCUT(fitu,4,400,cutXU,cutYU);
      cutU->SetName(Form("cutU%d",id));
      cutU->SetVarX("nebulaQPed");
      cutU->SetVarY("nebulaTOFUGP");
      cutU->Draw("SAME");

      c2->cd(2);
      chain->Draw(Form("nebulaTOFUGP:nebulaQPed>>huu%d(100,0,4000,500,-20,80)",id),Form("cutU%d&&nebulaID==%d&&abs(nebulaTOFUGP-%lf+8)<7",id,id,urange),"colz");
      huu = (TH2D*)gDirectory->Get(Form("huu%d",id));
      hpuu = (TProfile*)huu->ProfileX(Form("hpuu%d",id));
      hpuu->Fit("fitu","q","same");
      huu->GetYaxis()->SetRangeUser(urange-20,urange+5);
      huu->GetZaxis()->SetRangeUser(0.01,100);

      TCutG* cutUU = makeCUT(fitu,2,400,cutXU,cutYU);
      cutUU->SetName(Form("cutUU%d",id));
      cutUU->SetVarX("nebulaQPed");
      cutUU->SetVarY("nebulaTOFUGP");
      cutUU->Draw("SAME");

      c2->cd(3);
      chain->Draw(Form("nebulaTOFUGP:nebulaQPed>>huuu%d(100,0,4000,500,-20,80)",id),Form("cutUU%d&&nebulaID==%d&&abs(nebulaTOFUGP-%lf+8)<7",id,id,urange),"colz");
      huuu = (TH2D*)gDirectory->Get(Form("huuu%d",id));
      hpuuu = (TProfile*)huuu->ProfileX(Form("hpuuu%d",id));
      hpuuu->Fit("fitu","q","same");
      huuu->GetYaxis()->SetRangeUser(urange-20,urange+5);
      huuu->GetZaxis()->SetRangeUser(0.01,100);

      fout1<<fitu->GetParameter(0)<<" "<<fitu->GetParameter(1)<<endl;

      c2->cd(4);
      chain->Draw(Form("nebulaTOFDGP:nebulaQPed>>hd%d(100,0,4000,500,-20,80)",id),Form("nebulaID==%d&&abs(nebulaTOFDGP-%lf+8)<7",id,drange),"colz");
      hd = (TH2D*)gDirectory->Get(Form("hd%d",id));
      hpd = (TProfile*)hd->ProfileX(Form("hpd%d",id));
      TF1 *fitd = new TF1("fitd","[0]+[1]*TMath::Power(x,-0.5)",0,4000);
      fitd->SetParameter(1,40);
      fitd->SetParLimits(1,20,60);
      hpd->Fit("fitd","rq","same",10,500);

      hd->GetYaxis()->SetRangeUser(drange-20,drange+5);
      hd->GetZaxis()->SetRangeUser(0.01,100);

      //TCutG* makeCUT(TF1 *func, double width, int pointNum, double *cutX, double *cutY)

      double *cutXD;
      double *cutYD;
      TCutG* cutD = makeCUT(fitd,4,400,cutXD,cutYD);
      cutD->SetName(Form("cutD%d",id));
      cutD->SetVarX("nebulaQPed");
      cutD->SetVarY("nebulaTOFDGP");
      cutD->Draw("SAME");


      c2->cd(5);
      chain->Draw(Form("nebulaTOFDGP:nebulaQPed>>hdd%d(100,0,4000,500,-20,80)",id),Form("cutD%d&&nebulaID==%d&&abs(nebulaTOFDGP-%lf+8)<7",id,id,drange),"colz");
      hdd = (TH2D*)gDirectory->Get(Form("hdd%d",id));
      hpdd = (TProfile*)hdd->ProfileX(Form("hpdd%d",id));
      hpdd->Fit("fitd","q","same");
      hdd->GetYaxis()->SetRangeUser(drange-20,drange+5);
      hdd->GetZaxis()->SetRangeUser(0.01,100);

      TCutG* cutDD = makeCUT(fitd,2,400,cutXD,cutYD);
      cutDD->SetName(Form("cutDD%d",id));
      cutDD->SetVarX("nebulaQPed");
      cutDD->SetVarY("nebulaTOFDGP");
      cutDD->Draw("SAME");

      c2->cd(6);
      chain->Draw(Form("nebulaTOFDGP:nebulaQPed>>hddd%d(100,0,4000,500,-20,80)",id),Form("cutDD%d&&nebulaID==%d&&abs(nebulaTOFDGP-%lf+8)<7",id,id,drange),"colz");
      hddd = (TH2D*)gDirectory->Get(Form("hddd%d",id));
      hpddd = (TProfile*)hddd->ProfileX(Form("hpddd%d",id));
      hpddd->Fit("fitd","q","same");
      hddd->GetYaxis()->SetRangeUser(drange-20,drange+5);
      hddd->GetZaxis()->SetRangeUser(0.01,100);

      fout2<<fitd->GetParameter(0)<<" "<<fitd->GetParameter(1)<<endl;

      c1->Update();
      c2->Update();

      c1->Print(Form("./fig/NEBULA/slew_2nd_nebula_pol1_layer%d_1.pdf",layer));
      c2->Print(Form("./fig/NEBULA/slew_2nd_nebula_pol1_layer%d_2.pdf",layer));

    }

  c1->Print(Form("./fig/NEBULA/slew_2nd_nebula_pol1_layer%d_1.pdf]",layer));
  c2->Print(Form("./fig/NEBULA/slew_2nd_nebula_pol1_layer%d_2.pdf]",layer));
  fout1.close();
  fout2.close();
}






void save_hist()
{
  TFile *file = new TFile("./hist/NEBULA/slew_2nd_nebula.root","recreate");

  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  int nebulaNum;
  int* nebulaID = new int[144];
  double* nebulaTOF = new double[144];
  double* nebulaTOFUGP = new double[144];
  double* nebulaTOFDGP = new double[144];
  double* nebulaQPed = new double[144];
  double* nebulaQA = new double[144];

  chain->SetBranchAddress("nebulaNum",&nebulaNum);
  chain->SetBranchAddress("nebulaID",nebulaID);
  chain->SetBranchAddress("nebulaTOF",nebulaTOF);
  chain->SetBranchAddress("nebulaTOFUGP",nebulaTOFUGP);
  chain->SetBranchAddress("nebulaTOFDGP",nebulaTOFDGP);
  chain->SetBranchAddress("nebulaQPed",nebulaQPed);
  chain->SetBranchAddress("nebulaQA",nebulaQA);

  TH2 *h1;
  TH2 *h2;
  TH2 *h3;

  int id = 15;
  for(id = 1 ; id <= 120 ; id++)
    {
      cout<<"ID"<<id<<" is ongoing..."<<"\r";
      cout.flush();

      chain->Draw(Form("nebulaTOF:nebulaQPed>>h1id%d(4000,0,4000,500,-20,80)",id),Form("nebulaID==%d",id),"goff");
      chain->Draw(Form("nebulaTOFUGP:nebulaQPed>>h2id%d(4000,0,4000,500,-20,80)",id),Form("nebulaID==%d",id),"goff");
      chain->Draw(Form("nebulaTOFDGP:nebulaQPed>>h3id%d(4000,0,4000,500,-20,80)",id),Form("nebulaID==%d",id),"goff");

      h1 = (TH2D*)gDirectory->Get(Form("h1id%d",id));
      h2 = (TH2D*)gDirectory->Get(Form("h2id%d",id));
      h3 = (TH2D*)gDirectory->Get(Form("h3id%d",id));

      h1->Write();
      h2->Write();
      h3->Write();



    }
  file->Close();
}
