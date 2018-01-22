void slew_nebula()
{
  //22Ne beam
  //Al target
  //gamma ray production



  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  ofstream fout("./dat/nebula_slew.dat");
  
  TH2 *h2;
  TH1 *h3;
  TH2 *h4;
  TProfile *h5;
  TH2 *h6;
  TH2 *h7;
  TProfile *h8;
  TH2 *h9;
  
  TSpectrum *ts = new TSpectrum();
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->Print("./fig/slew_bad.pdf[");

  TCanvas *c2 = new TCanvas("c2","c2",1500,1000);
  c2->Divide(3,2);
  c2->Print("./fig/slew_cor.pdf[");

  TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
  c3->Divide(2,2);
  c3->Print("./fig/slew_cor_result.pdf[");
  
  Int_t numBAD=0;
  
  for(Int_t id=1;id<=120;id++)
    {
      c1->cd(1);
      chain->Draw(Form("nebulaTOF:nebulaQPed>>hi2d%d(1000,0,4000,200,-300,-200)",id),Form("nebulaID==%d",id),"colz");
      c1->cd(2);
      chain->Draw(Form("nebulaTOF>>hi%d(100,-300,-200)",id),Form("nebulaID==%d",id));

      h2 = (TH2D*)gDirectory->Get(Form("hi%d",id));
      h3 = (TH1D*)gDirectory->Get(Form("hi2d%d",id));
      h2->GetZaxis()->SetRangeUser(0.1,10);
      
      ts->Search(h2,1,"",0.02);

      Double_t center = ts->GetPositionX()[1];


            
      if(ts->GetNPeaks()!=2)
	{
	  numBAD++;
	  c1->Update();
	  c1->Print("./fig/slew_bad.pdf");
	  continue;
	}
      
     
	      
      c2->cd(1);
      chain->Draw(Form("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h4ID%d(200,0,0.4,400,-280,-230)",id),Form("nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,center),"colz");
      h4 = (TH2D*)gDirectory->Get(Form("h4ID%d",id));
      Double_t w4 = 5;
      Double_t d4 = h4->GetMean(2) - w4;
      Double_t u4 = h4->GetMean(2) + w4;
      h4->GetYaxis()->SetRangeUser(d4,u4);
      h4->GetZaxis()->SetRangeUser(0.1,10);
      
      c2->cd(2);
      chain->Draw(Form("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h5ID%d(100,0,0.4,400,-280,-230)",id),Form("nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,center),"prof");
      h5 = (TProfile*)gDirectory->Get(Form("h5ID%d",id));
      h5->GetYaxis()->SetRangeUser(d4,u4);
      TF1 *fit1 =new TF1("fit1","pol2",0.16,0.36);
      fit1->SetParameter(0,center);
      h5->Fit("fit1","rq","",0.16,0.36);


      //1st cut//
      Double_t cutX1[14]={0,};
      Double_t cutY1[14]={0,};

      Double_t c4 = h4->GetMean(1);
      Double_t s4 = h4->GetStdDev(1);

      for(Int_t i=0;i<7;i++)
	{
	  cutX1[i] = c4 - 3*s4 + i*s4;
	  cutY1[i] = fit1->Eval(cutX1[i]) + 0.5;

	  cutX1[i+7] = c4 + 3*s4 - i*s4;
	  cutY1[i+7] = fit1->Eval(cutX1[i+7]) - 1;
	}

      TCutG *cut1 = new TCutG(Form("cut1ID%d",id),14,cutX1,cutY1);
      c2->cd(1);
      cut1->SetVarX("TMath::Power(nebulaQPed,-0.25)");
      cut1->SetVarY("nebulaTOF");
      cut1->Draw("SAME");
      //1st cut//

      cout<<fit1->GetParameter(0)<<" "<<fit1->GetParameter(1)<<" "<<fit1->GetParameter(2)<<endl;
      Double_t slewPar0 = fit1->GetParameter(0);
      Double_t slewPar1 = fit1->GetParameter(1);
      Double_t slewPar2 = fit1->GetParameter(2);


      c2->cd(1);
      fit1->Draw("SAMEL");
      
      c2->cd(3);
      chain->Draw(Form("nebulaTOF-%lf*TMath::Power(nebulaQPed,-0.25)-%lf*TMath::Power(nebulaQPed,-0.5):TMath::Power(nebulaQPed,-0.25)>>h6ID%d(200,0,0.4,400,-280,-230)",slewPar1,slewPar2,id),Form("nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,center),"colz");
      h6 = (TH2D*)gDirectory->Get(Form("h6ID%d",id));
      Double_t w6 = w4;
      Double_t d6 = h6->GetMean(2) - w6;
      Double_t u6 = h6->GetMean(2) + w6;
      h6->GetYaxis()->SetRangeUser(d6,u6);
      h6->GetZaxis()->SetRangeUser(0.1,10);

      
      c2->cd(4);
      chain->Draw(Form("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h7ID%d(200,0,0.4,400,-280,-230)",id),Form("cut1ID%d&&nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,id,center),"goff");
      h7 = (TH2D*)gDirectory->Get(Form("h7ID%d",id));
      h7->GetYaxis()->SetRangeUser(d4,u4);
      h7->GetZaxis()->SetRangeUser(0.1,10);

      //2nd cut//
      c2->cd(5);
      chain->Draw(Form("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h8ID%d(100,0,0.4,400,-280,-230)",id),Form("cut1ID%d&&nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,id,center),"prof");
      h8 = (TProfile*)gDirectory->Get(Form("h8ID%d",id));
      h8->GetYaxis()->SetRangeUser(d4,u4);
      TF1 *fit2 = new TF1("fit2","pol2",0.16,0.36);
      fit2->SetParameters(fit1->GetParameters());
      h8->Fit("fit2","rq","",0.16,0.36);

      Double_t cutX2[14]={0,};
      Double_t cutY2[14]={0,};

      Double_t c7 = h7->GetMean(1);
      Double_t s7 = h7->GetStdDev(1);

      
      for(Int_t i=0;i<7;i++)
	{
	  cutX2[i] = c7 - 3*s7 + i*s7;
	  cutY2[i] = fit2->Eval(cutX2[i]) + 0.5;

	  cutX2[i+7] = c7 + 3*s7 - i*s7;
	  cutY2[i+7] = fit2->Eval(cutX2[i+7]) - 1;
	}

      TCutG *cut2 = new TCutG(Form("cut2ID%d",id),14,cutX2,cutY2);
      c2->cd(1);
      cut2->SetVarX("TMath::Power(nebulaQPed,-0.25)");
      cut2->SetVarY("nebulaTOF");
      cut2->SetLineColor(4);
      cut2->Draw("SAME");

      //2nd cut//
      
      Int_t nentry = chain->GetEntries(Form("cut2ID%d&&nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,id,center));
      if(nentry>100)
	cout<<"ID"<<id<<"'s slew correction is ongoing..."<<endl;
      else
	{
	  numBAD++;
	  c1->Update();
	  //c1->Print("./fig/slew_bad.pdf");
	  continue;
	}

      
      
      
      c2->cd(4);
      chain->Draw(Form("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h10ID%d(200,0,0.4,400,-280,-230)",id),Form("cut2ID%d&&nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,id,center),"colz");
      h7 = (TH2D*)gDirectory->Get(Form("h10ID%d",id));
      h7->GetYaxis()->SetRangeUser(d4,u4);
      h7->GetZaxis()->SetRangeUser(0.1,10);

      
      c2->cd(5);
      chain->Draw(Form("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h11ID%d(100,0,0.4,400,-280,-230)",id),Form("cut2ID%d&&nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,id,center),"prof");
      h8 = (TProfile*)gDirectory->Get(Form("h11ID%d",id));
      h8->GetYaxis()->SetRangeUser(d4,u4);
      Double_t fitl = h7->GetMean(1) - 2 * h7->GetStdDev(1);
      Double_t fitr = h7->GetMean(1) + 2 * h7->GetStdDev(1);
      TF1 *fit3 = new TF1("fit3","pol2",fitl,fitr);
      fit3->SetParameters(fit2->GetParameters());
      h8->Fit("fit3","rq","",fitl,fitr);
      
      cout<<fit3->GetParameter(0)<<" "<<fit3->GetParameter(1)<<" "<<fit3->GetParameter(2)<<endl;
      cout<<endl;
      slewPar0 = fit3->GetParameter(0);
      slewPar1 = fit3->GetParameter(1);
      slewPar2 = fit3->GetParameter(2);

      fout<<id<<" "<<slewPar0<<" "<<slewPar1<<" "<<slewPar2<<" "<<h4->GetEntries()<<endl;

      
      c2->cd(4);
      fit3->Draw("SAMEL");

      
      c2->cd(6);
      chain->Draw(Form("nebulaTOF-%lf*TMath::Power(nebulaQPed,-0.25)-%lf*TMath::Power(nebulaQPed,-0.5):TMath::Power(nebulaQPed,-0.25)>>h9ID%d(200,0,0.4,400,-280,-230)",slewPar1,slewPar2,id),Form("nebulaID==%d&&abs(nebulaTOF-%lf)<5",id,center),"colz");
      
      h9 = (TH2D*)gDirectory->Get(Form("h9ID%d",id));
      Double_t w9 = w4;
      Double_t d9 = h9->GetMean(2) - w9;
      Double_t u9 = h9->GetMean(2) + w9;
      h9->GetYaxis()->SetRangeUser(d9,u9);
      h9->GetZaxis()->SetRangeUser(0.1,10);

      c3->cd(1);
      h3->Draw("colz");
      c3->cd(2);
      h2->Draw();
      c3->cd(3);
      chain->Draw(Form("nebulaTOF-%lf*TMath::Power(nebulaQPed,-0.25)-%lf*TMath::Power(nebulaQPed,-0.5):nebulaQPed>>h12ID%d(1000,0,4000,200,-300,-200)",slewPar1,slewPar2,id),Form("nebulaID==%d",id),"colz");
      c3->cd(4);
      chain->Draw(Form("nebulaTOF-%lf*TMath::Power(nebulaQPed,-0.25)-%lf*TMath::Power(nebulaQPed,-0.5)>>h13ID%d(100,-300,-200)",slewPar1,slewPar2,id),Form("nebulaID==%d",id));
      
        
      c1->Update();
      c2->Update();
      c3->Update();
      c2->Print("./fig/slew_cor.pdf");
      c3->Print("./fig/slew_cor_result.pdf");
      // if(id==9)
      // 	break;
      
      //getchar();
      

    }
  c1->Print("./fig/slew_bad.pdf]");
  c2->Print("./fig/slew_cor.pdf]");
  c3->Print("./fig/slew_cor_result.pdf]");
  cout<<"Numbers of bad histogram : "<<numBAD<<endl;
      
}

void slew_additional(int i=0)
{
  /*
    Int_t numBad1[42] = {1,2,3,4,5,27,28,29,30
    ,31,32,33,34,35,36,57,58,59,60
    ,61,62,63,64,65,66,87,88,89,90
    ,91,92,93,94,95,96,97,100,102,106,118,119,120};

    //9 + 10 + 10 + 13 = 42
    */

  Int_t numBad1[37] = {1,2,3,4,5,27,28,29,30
		       ,31,32,33,34,35,36,57,58,59,60
		       ,61,62,63,64,65,66,87,88,89,90
		       ,95,96,97,100,102,106,118,119};
  //9 + 10 + 10 + 8 = 37

  
  // no gamma-ray line //
  //Int_t numBad2[] = {61,62,63,64,90,91,92,93,94,95,96,119,120};
  Int_t numBad2[5] = {91,92,93,94,120};
  // no gamma-ray line //


  
  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");


  TCanvas *c1 = new TCanvas(Form("c%d",i),"c1",1800,900);
  c1->Divide(3,2);

  TH2 *h0;
  TH2 *h1;
  TH1 *h2;
  TH2 *h3;
  TProfile *h4;
  TH1 *h5;

  c1->Print("./fig/slew_nebula_bad.pdf[");

  ofstream fout("./dat/slew_nebula_bad.dat");
  
  for(i=0;i<37;i++)
    {
     
      Int_t id = numBad1[i];
      Int_t layer;
      if(i<9)
	layer=1;
      else if(i<19)
	layer=2;
      else if(i<29)
	layer=3;
      else if(i<37)
	layer=4;
      else
	return;

      TFile *file = new TFile(Form("./cut/nebula_slew_add_layer%d.root",layer),"r");
      TCutG *cut = (TCutG*)file->Get(Form("cutID%d",id));

      cout<<"ID"<<id<<"'s slew correction is ongoing..."<<endl;

      c1->cd(1);
      chain->Draw("nebulaTOF:nebulaQPed>>h0(200,30,430,200,-270,-240)",Form("nebulaID==%d",id),"colz");
      h0 = (TH2D*)gDirectory->Get("h0");
      h0->GetZaxis()->SetRangeUser(0.1,10);
      cut->Draw("SAME"); 

  
      c1->cd(4);
      chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-270,-240)",Form("cutID%d&&nebulaID==%d",id,id),"colz");
      h1 = (TH2D*)gDirectory->Get("h1");
      h1->GetZaxis()->SetRangeUser(0.1,10);
 
  
      c1->cd(3);
      chain->Draw("nebulaTOF>>h2(100,-270,-240)",Form("cutID%d&&nebulaID==%d",id,id));
      h2 = (TH1D*)gDirectory->Get("h2");

  
      c1->cd(2);
      chain->Draw("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h3(200,0.1,0.5,200,-270,-240)",Form("cutID%d&&nebulaID==%d",id,id),"colz");
      h3 = (TH2D*)gDirectory->Get("h3");
      h3->GetZaxis()->SetRangeUser(0.1,10);

      TH1 *h33 = (TH1D*)h3->ProjectionX("h33");
  
      Double_t left  = h33->GetBinCenter(h33->FindFirstBinAbove(0,1));
      Double_t right = h33->GetBinCenter(h33->FindLastBinAbove(0,1));
  
  
      c1->cd(5);
      chain->Draw("nebulaTOF:TMath::Power(nebulaQPed,-0.25)>>h4(100,0.1,0.5,200,-270,-240)",Form("cutID%d&&nebulaID==%d",id,id),"prof");
      h4 = (TProfile*)gDirectory->Get("h4");
      h4->GetYaxis()->SetRangeUser(-270,-240);
      TF1 *h44 = new TF1("h44","pol2");
      h4->Fit("h44","rq","",left,right);

      TF1 *h444 = new TF1("h444","[0] + [1]*TMath::Power(x,-0.25) + [2]*TMath::Power(x,-0.5)",40,400);
      h444->SetParameters(h44->GetParameters());
      h444->SetLineColor(2);

      Double_t par0 = h44->GetParameter(0);
      Double_t par1 = h44->GetParameter(1);
      Double_t par2 = h44->GetParameter(2);
  
      c1->cd(6);
      Double_t l5 = h44->GetParameter(0)-15;
      Double_t r5 = h44->GetParameter(0)+15;
      chain->Draw(Form("nebulaTOF-%lf*TMath::Power(nebulaQPed,-0.25)-%lf*TMath::Power(nebulaQPed,-0.5)>>h5(100,%lf,%lf)",par1,par2,l5,r5),Form("cutID%d&&nebulaID==%d",id,id));
      h5 = (TH1D*)gDirectory->Get("h5");

      if(par2>0)
	fout<<id<<" "<<par0<<" "<<par1<<" "<<par2<<" "<<h1->GetEntries()<<endl;  
  
      c1->cd(2);
      h44->Draw("SAME");
      c1->cd(1);
      h444->Draw("SAME");
      c1->cd(4);
      h444->Draw("SAME");
      c1->Update();

      c1->Print("./fig/slew_nebula_bad.pdf");
    }
  fout.close();
  c1->Print("./fig/slew_nebula_bad.pdf]");
}

void slew_add_rep(int layer=1)
{
  //9 + 10 + 10 + 8 = 37
  int Nmin = 0;
  int Nmax = 9;
  if(layer==1)
    {
      Nmin = 0;
      Nmax = 9;
    }
  else if(layer==2)
    {
      Nmin = 9;
      Nmax = 19;
    }
  else if(layer==3)
    {
      Nmin = 19;
      Nmax = 29;
    }
  else if(layer==4)
    {
      Nmin = 29;
      Nmax = 37;
    }
  for(Int_t i=Nmin;i<Nmax;i++)
    {
      slew_additional(i);
      getchar();
    }
}

void draw_additional(int i=0)
{
  /*
    Int_t numBad1[7] = {1,2,3,4,5,29,30};
    Int_t numBad2[9] = {31,32,33,34,35,36,58,59,60};
    Int_t numBad3[9] = {61,62,63,64,65,66,88,89,90};
    Int_t numBad4[9] = {91,92,93,94,95,96,118,119,120};
    Int_t numBad5[9] = {27,28,57,87,97,99,100,102,106};
  */

  /*
    Int_t numBad1[42] = {1,2,3,4,5,27,28,29,30
    ,31,32,33,34,35,36,57,58,59,60
    ,61,62,63,64,65,66,87,88,89,90
    ,91,92,93,94,95,96,97,100,102,106,118,119,120};
    //9 + 10 + 10 + 13 = 42
    */

  Int_t numBad1[37] = {1,2,3,4,5,27,28,29,30
		       ,31,32,33,34,35,36,57,58,59,60
		       ,61,62,63,64,65,66,87,88,89,90
		       ,95,96,97,100,102,106,118,119};
  //9 + 10 + 10 + 8 = 37

  
  // no gamma-ray line //
  //Int_t numBad2[] = {61,62,63,64,90,91,92,93,94,95,96,119,120};
  Int_t numBad2[5] = {91,92,93,94,120};
  // no gamma-ray line //

  
  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");


  Int_t id = numBad1[i];


  chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-280,-220)",Form("nebulaID==%d",id),"colz");

  
}

void slew_nebula_old()
{
  //22Ne beam
  //Al target
  //gamma ray production



  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  ofstream fout("./dat/nebula_slew.dat");
  
  TH2 *h2;
  TH1 *h3;
  TH2 *h4;
  TProfile *h5;
  TH2 *h6;
  TSpectrum *ts = new TSpectrum();
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  TCanvas *c2 = new TCanvas("c2","c2",1800,600);
  c1->Divide(2,1);
  c1->Print("./fig/slew_old_bad.pdf[");
  c2->Divide(3,1);
  c2->Print("./fig/slew_old_cor.pdf[");

  Int_t numBAD=0;
  
  for(Int_t id=1;id<=120;id++)
    {
      //cout<<"ID : "<<id<<endl;
      c1->cd(1);
      chain->Draw(Form("nebulaTA:nebulaQPed>>hi2d%d(400,0,300,200,5,100)",id),Form("nebulaID==%d",id),"colz");
      //chain->Draw(Form("nebulaTU:nebulaQUPed>>hi2d%d(400,0,300,200,5,100)",id),Form("nebulaID==%d",id),"colz");
      c1->cd(2);
      chain->Draw(Form("nebulaTA>>hi%d(100,5,100)",id),Form("nebulaID==%d",id),"colz");
      //chain->Draw(Form("nebulaTU>>hi%d(100,5,100)",id),Form("nebulaID==%d",id),"colz");

      h2 = (TH2D*)gDirectory->Get(Form("hi%d",id));
      h3 = (TH1D*)gDirectory->Get(Form("hi2d%d",id));
      h2->GetZaxis()->SetRangeUser(0.1,10);
      //h3->GetZaxis()->SetRangeUser(0.1,10);
      
      ts->Search(h2,1,"",0.02);

      /*cout<<"peak numbers : "<<ts->GetNPeaks()<<endl;
	cout<<ts->GetPositionX()[0]<<" "<<ts->GetPositionX()[1]<<endl;
      */
      if(ts->GetNPeaks() == 2)
	cout<<"ID"<<id<<"'s slew correction is ongoing..."<<endl;
      else
	{
	  numBAD++;
	  c1->Update();
	  c1->Print("./fig/slew_old_bad.pdf");
	  continue;
	}      
      c2->cd(1);
      Double_t center = ts->GetPositionX()[1];
      chain->Draw(Form("nebulaTA:nebulaQPed>>h4ID%d(200,0,300,200,5,100)",id),Form("nebulaID==%d&&abs(nebulaTA-%lf-2)<3",id,center),"colz");
      h4 = (TH2D*)gDirectory->Get(Form("h4ID%d",id));
      h4->GetYaxis()->SetRangeUser(30,70);
      h4->GetZaxis()->SetRangeUser(0.1,10);

      c2->cd(2);
      chain->Draw(Form("nebulaTA:nebulaQPed>>h5ID%d(50,0,300,200,5,100)",id),Form("nebulaID==%d&&abs(nebulaTA-%lf-2)<3",id,center),"prof");
      h5 = (TProfile*)gDirectory->Get(Form("h5ID%d",id));
      h5->GetYaxis()->SetRangeUser(30,70);
      TF1 *fit = new TF1("fit","[0]+[1]*1./TMath::Power(x,0.5)",50,250);
      fit->SetParameter(0,center);
      fit->SetParameter(1,40);
      //fit->SetParLimits(1,20,80);
      h5->Fit("fit","rq","",50,250);

      
      cout<<fit->GetParameter(0)<<" "<<fit->GetParameter(1)<<endl;
      cout<<endl;
      Double_t slewPar1 = fit->GetParameter(1);

      fout<<id<<" "<<slewPar1<<" "<<h4->GetEntries()<<endl;
      
      c2->cd(1);
      fit->Draw("SAMEL");


      c2->cd(3);
      chain->Draw(Form("nebulaTA-%lf/sqrt(nebulaQPed):nebulaQPed>>h6ID%d(200,0,300,200,5,100)",slewPar1,id),Form("nebulaID==%d&&abs(nebulaTA-%lf-2)<3",id,center),"colz");
      //chain->Draw(Form("nebulaTU-%lf/sqrt(nebulaQUPed):nebulaQUPed>>h6ID%d(200,0,300,200,5,100)",slewPar1,id),Form("nebulaID==%d&&abs(nebulaTU-%lf-2)<3",id,center),"colz");
      h6 = (TH2D*)gDirectory->Get(Form("h6ID%d",id));
      h6->GetYaxis()->SetRangeUser(30,70);
      h6->GetZaxis()->SetRangeUser(0.1,10);
      
      c1->Update();
      c2->Update();
      c2->Print("./fig/slew_old_cor.pdf");
      // if(id==9)
      // 	break;
      
      //getchar();
      
      
    }
  c1->Print("./fig/slew_old_bad.pdf]");
  c2->Print("./fig/slew_old_cor.pdf]");
  cout<<"Numbers of bad histogram : "<<numBAD<<endl;
}




void save_hist()
{

  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");


  TH2 *h1;
  chain->Draw("nebulaTA:nebulaID>>h1(144,0.5,144.5,1000,-50,250)","","goff");
  h1 = (TH2D*)gDirectory->Get("h1");

  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);

  TH2 *h2;
  TH2 *h3;
  TSpectrum *ts = new TSpectrum();

  c1->Print("./fig/slew.pdf[");
  
  for(Int_t id=1;id<=120;id++)
    {
      c1->cd(1);
      chain->Draw(Form("nebulaTA:nebulaQPed>>hi2d%d(400,0,300,200,5,100)",id),Form("nebulaID==%d",id),"colz");
      c1->cd(2);
      chain->Draw(Form("nebulaTA>>hi%d(100,5,100)",id),Form("nebulaID==%d",id),"colz");

      h2 = (TH2D*)gDirectory->Get(Form("hi%d",id));
      h3 = (TH2D*)gDirectory->Get(Form("hi2d%d",id));
      //h2->GetZaxis()->SetRangeUser(0.1,10);
      h3->GetZaxis()->SetRangeUser(0.1,10);
      
      ts->Search(h2);
      c1->Update();
      c1->Print("./fig/slew.pdf");
      //getchar();
      
      
    }
  c1->Print("./fig/slew.pdf]");
  
}
void check_par()
{
  ifstream fin1("./dat/nebula_slew.dat");
  ifstream fin2("./dat/slew_nebula_bad.dat");
  ofstream fout("./dat/nebula_slew_new.dat");
  double id[200]={0,};
  double par0[200]={0,};
  double par1[200]={0,};
  double par2[200]={0,};
  double eventNum[200]={0,};
  int index[200]={0,};
  int i=0;

  TTree *tree = new TTree("t1","t1");
  double ID;
  double p1;
  double p2;
  double entry;
  /*
    tree->Branch("ID",&ID);
    tree->Branch("p1",&p1);
    tree->Branch("p2",&p2);
    tree->Branch("entry",&entry);
  */

  tree->Branch("ID",&ID,"ID/D");
  tree->Branch("p1",&p1,"p1/D");
  tree->Branch("p2",&p2,"p2/D");
  tree->Branch("entry",&entry,"entry/D");
  
  while(1)
    {
      if(fin1.fail())
	{
	  i--;
	  break;
	}
      fin1>>id[i]>>par0[i]>>par1[i]>>par2[i]>>eventNum[i];
      ID=id[i];
      p1=par1[i];
      p2=par2[i];
      entry=eventNum[i];

      tree->Fill();
      //cout<<i<<" "<<id[i]<<" "<<par1[i]<<" "<<par2[i]<<" "<<eventNum[i]<<endl;
      i++;
    }

  double id_last[120];
  double par0_last[120];
  double par1_last[120];
  double par2_last[120];
  double eventNum_last[120];



  while(1)
    {
      bool flag = false;
      if(fin2.fail())
	{
	  i--;
	  break;
	}
      double dummyID, dummyP0, dummyP1, dummyP2, dummyEN;
      fin2>>dummyID>>dummyP0>>dummyP1>>dummyP2>>dummyEN;
      for(int j=0;j<200;j++)
	{
	  if(dummyID == id[j])
	    {
	      par0[j] = dummyP0;
	      par1[j] = dummyP1;
	      par2[j] = dummyP2;
	      eventNum[j] = dummyEN;
	      flag = true;
	      break;
	    }
	}
      if(flag == true)
	continue;
      id[i] = dummyID;
      par0[i] = dummyP0;
      par1[i] = dummyP1;
      par2[i] = dummyP2;
      eventNum[i] = dummyEN;
      dummyID=0;dummyP0=0;dummyP1=0;dummyP2=0;dummyEN=0;
      //cout<<i<<" "<<id[i]<<" "<<par1[i]<<" "<<par2[i]<<endl;
      i++;
    }

  int NN = 0;
  double avg_par0 = 0;
  double avg_par1 = 0;
  double avg_par2 = 0;
  TMath::Sort(200,id,index);

  int event_cri = 500;
  
  for(i=0;i<120;i++)
    {
      id_last[i] = id[index[i]];
      par0_last[i] = par0[index[i]];
      par1_last[i] = par1[index[i]];
      par2_last[i] = par2[index[i]];
      eventNum_last[i] = eventNum[index[i]];

      if(id_last[i] == 0)
	id_last[i] = 9999;
      
      if(eventNum_last[i]>event_cri)
	{
	  NN += eventNum[i];
	  avg_par0 += par0_last[i] * eventNum[i] ;
	  avg_par1 += par1_last[i] * eventNum[i] ;
	  avg_par2 += par2_last[i] * eventNum[i] ;
	  //cout<<i<<" "<<id_last[i]<<" "<<par1_last[i]<<" "<<par2_last[i]<<" "<<eventNum_last[i]<<endl;
	}
      //cout<<i<<" "<<id_last[i]<<" "<<par1_last[i]<<" "<<par2_last[i]<<" "<<eventNum_last[i]<<endl;
    }

  avg_par0 = avg_par0/NN;
  avg_par1 = avg_par1/NN;
  avg_par2 = avg_par2/NN;

  cout<<endl;
  cout<<"# of entry>"<<event_cri<<" is "<<NN<<". "<<endl;
  cout<<"Average of p0 : "<<avg_par0<<endl;
  cout<<"Average of p1 : "<<avg_par1<<endl;
  cout<<"Average of p2 : "<<avg_par2<<endl;

  int index2[120];
  TMath::Sort(120,id_last,index2,false);
  for(i=0;i<120;i++)
    {
      id[i] = id_last[index2[i]];
      par0[i] = par0_last[index2[i]];
      par1[i] = par1_last[index2[i]];
      par2[i] = par2_last[index2[i]];
      eventNum[i] = eventNum_last[index2[i]];
      if(id[i] == 9999)
	{
	  par0[i] = avg_par0;
	  par1[i] = avg_par1;
	  par2[i] = avg_par2;
	}
      //cout<<i+1<<" "<<id[i]<<" "<<par1[i]<<" "<<par2[i]<<" "<<eventNum[i]<<endl;
    }

  int endIndex;
  std::vector<double> id_add;
  for(i=0;i<119;i++)
    {
      int diff = id[i+1] - id[i];
      if(diff == 1)
	continue;
      else if(diff > 100)
	endIndex = i;
      else if(id[i+1] == id[i])
	continue;
      else
	{
	  for(int j=1;j<diff;j++)
	    {
	      id_add.push_back(id[i]+j);
	      //cout<<id[i]+j<<endl;
	    }
	}
    }

  if(id[119] != 120)
    id[119] = 120;
   
  for(i=0;i<id_add.size();i++)
    {
      id[i+endIndex+1] = id_add.at(i);
      //cout<<id[i+endIndex+1]<<endl;
    }

  TMath::Sort(120,id,index2,false);
  for(i=0;i<120;i++)
    {
      id_last[i] = id[index2[i]];
      par0_last[i] = par0[index2[i]];
      par1_last[i] = par1[index2[i]];
      par2_last[i] = par2[index2[i]];
      eventNum_last[i] = eventNum[index2[i]];
      cout<<i+1<<" "<<id_last[i]<<" "<<par0_last[i]<<" "<<par1_last[i]<<" "<<par2_last[i]<<" "<<eventNum_last[i]<<endl;
      fout<<id_last[i]<<" "<<par0_last[i]<<" "<<par1_last[i]<<" "<<par2_last[i]<<" "<<eventNum_last[i]<<endl;
    }
  
  for(i=0;i<120;i++)
    {
      id[i] = id_last[i];
      par0[i] = par0_last[i];
      par1[i] = par1_last[i];
      par2[i] = par2_last[i];
      eventNum[i] = eventNum_last[i];
    }


  double u_range1 = 50;
  double d_range1 = -300;

  double u_range2 = 500;
  double d_range2 = -100;
  

  TGraph *gr1 = new TGraph(i,id,par1);
  gr1->SetMarkerStyle(24);
  gr1->SetTitle("Par1 VS ID;ID");
  gr1->GetXaxis()->SetRangeUser(0,120);
  gr1->GetYaxis()->SetRangeUser(d_range1,u_range1);
  
  TGraph *gr2 = new TGraph(i,id,par2);
  gr2->SetMarkerStyle(24);
  gr2->SetTitle("Par2 VS ID;ID");
  gr2->GetXaxis()->SetRangeUser(0,120);
  gr2->GetYaxis()->SetRangeUser(d_range2,u_range2);

  /*
    TGraph *gr3 = new TGraph(i,eventNum,par1);
    gr3->SetMarkerStyle(20);
    gr3->SetTitle("Par1 VS Entry number;Entry number");
  
    TGraph *gr4 = new TGraph(i,eventNum,par2);
    gr4->SetMarkerStyle(20);
    gr4->SetTitle("Par2 VS Entry number;Entry number");
  */

  
  TLine *l11 = new TLine(30,d_range1,30,u_range1);
  l11->SetLineStyle(2);
  l11->SetLineWidth(2);
  TLine *l12 = new TLine(60,d_range1,60,u_range1);
  l12->SetLineStyle(2);
  l12->SetLineWidth(2);
  TLine *l13 = new TLine(90,d_range1,90,u_range1);
  l13->SetLineStyle(2);
  l13->SetLineWidth(2);

  TLine *l14 = new TLine(0,avg_par1,120,avg_par1);
  l14->SetLineStyle(2);
  l14->SetLineWidth(2);
  l14->SetLineColor(2);

  TLine *l21 = new TLine(30,d_range2,30,u_range2);
  l21->SetLineStyle(2);
  l21->SetLineWidth(2);
  TLine *l22 = new TLine(60,d_range2,60,u_range2);
  l22->SetLineStyle(2);
  l22->SetLineWidth(2);
  TLine *l23 = new TLine(90,d_range2,90,u_range2);
  l23->SetLineStyle(2);
  l23->SetLineWidth(2);

  TLine *l24 = new TLine(0,avg_par2,120,avg_par2);
  l24->SetLineStyle(2);
  l24->SetLineWidth(2);
  l24->SetLineColor(2);


  
  TCanvas *c1 = new TCanvas("c1","c1",1600,600);
  c1->Divide(2,1);
  c1->cd(1);
  gr1->Draw("AP");
  l11->Draw("SAME");
  l12->Draw("SAME");  
  l13->Draw("SAME");
  l14->Draw("SAME");
  
  c1->cd(2);
  gr2->Draw("AP");
  l21->Draw("SAME");
  l22->Draw("SAME");  
  l23->Draw("SAME");
  l24->Draw("SAME");
  
  c1->Print("./fig/nebula_slew_par_check.pdf");
  
  fin1.close();
  fin2.close();
  fout.close();
}



void compare_res(int id)
{
  std::ofstream fout("./dat/nebula_slew_resolution.dat");
  fout<<setw(4)<<"id"<<"\t"<<setw(8)<<"entry"<<"\t"<<setw(8)<<"center"<<"\t"<<setw(8)<<"resolution"<<endl;
  std::ifstream fin2("./dat/nebula_slew_new.dat");
  Double_t nebula_slew_par0[120]={0,};
  Double_t nebula_slew_par1[120]={0,};
  Double_t nebula_slew_par2[120]={0,};
  Int_t nebula_gamma_entry[120]={0,};
  Double_t dummy;
  for(Int_t i=0;i<120;i++)
    fin2>>dummy>>nebula_slew_par0[i]>>nebula_slew_par1[i]>>nebula_slew_par2[i]>>nebula_gamma_entry[i];
  fin2.close();

  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->Divide(2,2);

  c1->Print("./fig/slew_nebula_resolution.pdf[");
  
  for(id = 0;id<120;id++)
    {
      //int id = 14;
      if(nebula_gamma_entry[id]<500)
	continue;

      //double offset = nebula_slew_par0[id];
      double offset = 0;
      TF1 *f1 = new TF1("f1","pol2");
      f1->SetParameter(0,nebula_slew_par0[id]);
      f1->SetParameter(1,nebula_slew_par1[id]);
      f1->SetParameter(2,nebula_slew_par2[id]);

      double urange1 = nebula_slew_par0[id]+8;
      double drange1 = nebula_slew_par0[id]-4;
      double urange2 = 8;
      double drange2 = -4;

      
      c1->cd(1);
      chain->Draw(Form("nebulaTOF:nebulaQPed>>h1(1000,0,2000,100,%lf,%lf)",drange1,urange1),Form("nebulaID==%d",id+1),"colz");
      TH2 *h1 = (TH2D*)gDirectory->Get("h1");
      h1->GetZaxis()->SetRangeUser(0.1,10);

      c1->cd(2);
      chain->Draw(Form("nebulaTOF>>h2(100,%lf,%lf)",drange1,urange1),Form("nebulaID==%d",id+1));
      TH1 *h2 = (TH1D*)gDirectory->Get("h2");

  
      c1->cd(3);
      chain->Draw(Form("nebulaTOFSlew:nebulaQPed>>h3(1000,0,2000,100,%lf,%lf)",drange2,urange2),Form("nebulaID==%d",id+1),"colz");
      TH2 *h3 = (TH2D*)gDirectory->Get("h3");
      h3->GetZaxis()->SetRangeUser(0.1,10);

      c1->cd(4);
      chain->Draw(Form("nebulaTOFSlew>>h4(100,%lf,%lf)",drange2,urange2),Form("nebulaID==%d",id+1));
      TH1 *h4 = (TH1D*)gDirectory->Get("h4");
      //h4->Fit("gaus","q0");
      TF1 *fit1 = new TF1("fit1","gaus",drange2,urange2);
      fit1->SetParameter(0,100);
      fit1->SetParameter(1,offset);
      fit1->SetParameter(2,0.3);
      fit1->SetParLimits(2,0,0.5);
      fit1->SetLineColor(1);
      fit1->SetLineStyle(2);
      h4->Fit("fit1","rq0","",-2,2);
      fit1 = (TF1*)h4->GetFunction("fit1");
      double lrange = fit1->GetParameter(1) - 2 * fit1->GetParameter(2);
      double rrange = fit1->GetParameter(1) + 2 * fit1->GetParameter(2);
      TF1 *fit2 = new TF1("fit2","gaus",lrange,rrange);
      fit2->SetParameters(fit1->GetParameters());
      h4->Fit("fit2","rq","",lrange,rrange);
      fout<<setw(4)<<id+1<<"\t"<<setw(8)<<nebula_gamma_entry[id]<<"\t"<<setw(8)<<fit2->GetParameter(1)<<"\t"<<setw(8)<<fit2->GetParameter(2)<<endl;

      c1->Update();
      c1->Print("./fig/slew_nebula_resolution.pdf");
      //getchar();
    }
  c1->Print("./fig/slew_nebula_resolution.pdf]");
}
