// Fitting Function (Klein - Nishina) //
double kn(double *x, double *par)
{
  double Eg   = 4.44;  // in MeV
  double C    = par[0];// total 
  double gain = par[1];// in MeV/Ch

  double xx = x[0] * gain / Eg;
  double AoB = 0.511 / Eg;

  double temp = 1 - AoB * xx/(1-xx);

  double result;

  double compEdgeRatio = 2 / (AoB + 2);

  if (xx < compEdgeRatio)
    result = C * (-xx + 1/(1-xx) + temp*temp);
  else
    result = 0;

  return result; 
}

double knexp(double *x, double *par)
{
  double Eg   = 4.44;  // in MeV
  double C    = par[0];// total 
  double gain = par[1];// in MeV/Ch

  double xx = x[0] * gain / Eg;
  double AoB = 0.511 / Eg;

  double temp = 1 - AoB * xx/(1-xx);

  double result;

  double compEdgeRatio = 2 / (AoB + 2);

  if (xx < compEdgeRatio)
    result = C * (-xx + 1/(1-xx) + temp*temp);
  else
    result = 0;

  // exponential back ground
  result  = result + TMath::Exp( par[2] + par[3] * xx );

  return result;
  
}


double knexp_gaus(double *x, double *par)
{
  double xx;
  double fkn;

  // convolution control
  int np = 1000; // # of steps
  double sc = 5; // gaussian sigmas +-

  // convolution integral
  double xlow = x[0] - sc * par[4];
  double xup  = x[0] + sc * par[4];

  double step = (xup - xlow) / np;

  double sum = 0;

  for(int i=0 ; i< np ; i++)
    {
      xx = xlow + (i - 0.5) * step;
      fkn = knexp(&xx,par);
      sum = sum + fkn * TMath::Gaus(x[0],xx,par[4],kTRUE);
    }

  return (sum * step);
}

void AmBePeak()
{
  //TFile *file = new TFile(Form("./root/nebula%04d.root.NEBULA",runNum));
  //TTree *tree = (TTree*)file->Get("NEBULA");

  TFile *file1 = new TFile("./root/nebula0030.root.NEBULA");
  TTree *tree1 = (TTree*)file1->Get("NEBULA");
  TFile *file2 = new TFile("./root/nebula0031.root.NEBULA");
  TTree *tree2 = (TTree*)file2->Get("NEBULA");
  TFile *file3 = new TFile("./root/nebula0032.root.NEBULA");
  TTree *tree3 = (TTree*)file3->Get("NEBULA");

  TFile *file4 = new TFile("./root/nebula0033.root.NEBULA");
  TTree *tree4 = (TTree*)file4->Get("NEBULA");
  TFile *file5 = new TFile("./root/nebula0034.root.NEBULA");
  TTree *tree5 = (TTree*)file5->Get("NEBULA");
  TFile *file6 = new TFile("./root/nebula0035.root.NEBULA");
  TTree *tree6 = (TTree*)file6->Get("NEBULA");

  //TChain *chain = new TChain("NEBULA","NEBULA");
  //chain->Add("./root/nebula003[0-5].root.NEBULA");
  
  TCanvas *c1 = new TCanvas("c1","c1",1800,400);
  c1->Divide(3,1);

  TH1 *h1;
  TH1 *h2;
  TH1 *h3;
  TH1 *h4;
  TH1 *h5;
  TH1 *h6;

  int number[6]={0,};
  int index[6]={0,};
  TString str[6];

  str[0] = "front right";
  str[1] = "front center";
  str[2] = "front left";	

  str[3] = "rear right";
  str[4] = "rear center";
  str[5] = "rear left";	

  //ID 1~10  : nebula0032
  //ID 11~18 : nebula0031
  //ID 19~30 : nebula0030

  //ID 31~42 : nebula0032
  //ID 43~48 : nebula0031
  //ID 49~60 : nebula0030
  
  //ID 61~69 : nebula0033
  //ID 70~90 : nebula0035
  //ID 91~100 : nebula0033
  //ID 101~120 : nebula0035
  
  int runNum[6] = {30,31,32,33,34,35}; 

  for(int id = 1 ; id<=120 ; id+=1)
    {
      c1->cd(1);
      tree1->Draw(Form("nebulaQUPed>>h1%d(400,0,400)",id),Form("nebulaID==%d",id),"goff");
      h1 = (TH1D*)gDirectory->Get(Form("h1%d",id));
      c1->cd(2);
      tree2->Draw(Form("nebulaQUPed>>h2%d(400,0,400)",id),Form("nebulaID==%d",id),"goff");
      h2 = (TH1D*)gDirectory->Get(Form("h2%d",id));
      c1->cd(3);
      tree3->Draw(Form("nebulaQUPed>>h3%d(400,0,400)",id),Form("nebulaID==%d",id),"goff");
      h3 = (TH1D*)gDirectory->Get(Form("h3%d",id));

      c1->cd(1);
      tree4->Draw(Form("nebulaQUPed>>h4%d(400,0,400)",id),Form("nebulaID==%d",id),"goff");
      h4 = (TH1D*)gDirectory->Get(Form("h4%d",id));
      c1->cd(2);
      tree5->Draw(Form("nebulaQUPed>>h5%d(400,0,400)",id),Form("nebulaID==%d",id),"goff");
      h5 = (TH1D*)gDirectory->Get(Form("h5%d",id));
      c1->cd(3);
      tree6->Draw(Form("nebulaQUPed>>h6%d(400,0,400)",id),Form("nebulaID==%d",id),"goff");
      h6 = (TH1D*)gDirectory->Get(Form("h6%d",id));

      number[0] = h1->GetEntries();
      number[1] = h2->GetEntries();
      number[2] = h3->GetEntries();
      number[3] = h4->GetEntries();
      number[4] = h5->GetEntries();
      number[5] = h6->GetEntries();
      
      c1->Update();
      cout<<"ID"<<id<<" ";
      for(int i=0;i<=5;i++)
	cout<<number[i]<<" ";
      cout<<endl;
      cout<<"ID"<<id<<" ";
      TMath::Sort(6,number,index);
      cout<<number[index[0]]<<endl;
      cout<<"Run #: ";
      cout<<runNum[index[0]]<<endl;
      cout<<"Position : "<<str[index[0]]<<endl;
      cout<<endl;
      
      if(id == 30 || id == 60 || id == 90 || id == 120)
	cout<<endl;
	 
    }

}




  void make2Dqdc()
  {
    TChain *chain = new TChain("NEBULA","NEBULA");
    chain->Add("./root/nebula003[6-8].root.NEBULA");
    TH2 *hist_qu = new TH2D("hist_qu","hist_qu",120,0.5,120.5,2000,0,4000);
    TH2 *hist_qd = new TH2D("hist_qd","hist_qd",120,0.5,120.5,2000,0,4000);

    chain->Draw("nebulaQUPed:nebulaID>>hist_qu","","goff");
    hist_qu = (TH2D*)gDirectory->Get("hist_qu");
    chain->Draw("nebulaQDPed:nebulaID>>hist_qd","","goff");
    hist_qd = (TH2D*)gDirectory->Get("hist_qd");

    TFile *file = new TFile("./hist/qdc_nebula.root","recreate");
    hist_qu->Write();
    hist_qd->Write();
    file->Close();
  
  }


  void findCosmicPeak()
  {

    TFile *file = new TFile("./hist/qdc_nebula.root","r");

    TH2 *hist_qu = (TH2D*)file->Get("hist_qu");
    TH2 *hist_qd = (TH2D*)file->Get("hist_qd");

  
    TH1 *h1;
    TH1 *h2;
    TF1 *f1;
    TF1 *f2;
  
    TF1 *cos_peak = new TF1("cos_peak","landau(0)+expo(3)");
    TF1 *cos_landau = new TF1("cos_landau","landau",0,4000);

    double fit_down = 400;
    double fit_up   = 2000; 

    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1);

    ofstream fout1("./dat/nebula_cosmic_peak_u.dat");
    ofstream fout2("./dat/nebula_cosmic_peak_d.dat");

    c1->Print("./fig/qdc_nebula.pdf[");
  
    //int id = 1;

    for(int id = 1; id<=120; id++)
      {
  
	c1->cd(1);
	h1 = hist_qu->ProjectionY(Form("qu%d",id),id,id);
	h1->SetTitle(Form("QDC ID%d Up;Q_{Ped} (channel)",id));
	h1->GetXaxis()->SetRangeUser(0,fit_up);
	h1->Fit("landau","q","",fit_down,fit_up);
	f1 = h1->GetFunction("landau");
	cos_peak->SetParameters(f1->GetParameters());
	cos_peak->SetParameter(3,5);
	cos_peak->SetParameter(4,-2E-2);
	h1->Fit(cos_peak,"q","",fit_down,fit_up);
	fout1<<cos_peak->GetParameter(1)<<" "<<cos_peak->GetParameter(2)<<endl;
  
	c1->cd(2);
	h2 = hist_qd->ProjectionY(Form("qd%d",id),id,id);
	h2->SetTitle(Form("QDC ID%d Down;Q_{Ped} (channel)",id));
	h2->GetXaxis()->SetRangeUser(0,fit_up);
	h2->Fit("landau","q","",fit_down,fit_up);
	f2 = h2->GetFunction("landau");
	cos_peak->SetParameters(f2->GetParameters());
	cos_peak->SetParameter(3,5);
	cos_peak->SetParameter(4,-2E-2);
	h2->Fit(cos_peak,"q","",fit_down,fit_up);
	fout2<<cos_peak->GetParameter(1)<<" "<<cos_peak->GetParameter(2)<<endl;
  
	c1->Print("./fig/qdc_nebula.pdf");

	cout<<"ID"<<id<<"is ongoing..."<<endl;

      }

    c1->Print("./fig/qdc_nebula.pdf]");
    fout1.close();
    fout2.close();
 
  }
