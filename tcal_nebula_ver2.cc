void find_peak()
{
  TFile *file =new TFile("./root/run0303.root.NEBULA");
  TTree *tree = (TTree*)file->Get("NEBULA");

  TH1 *ht1 = new TH1D("ht1","ht1",4500,0.5,4500.5);
  TH1 *ht2 = new TH1D("ht2","ht2",4500,0.5,4500.5);
  TSpectrum *ts1 = new TSpectrum();
  TSpectrum *ts2 = new TSpectrum();

  
  
  TH1 *hi1;
  for(Int_t id=1;id<=144;id++)
    {
      ofstream fout1(Form("./dat/NEBULA/time_vs_channel_u_id%03d.dat",id));
      ofstream fout2(Form("./dat/NEBULA/time_vs_channel_d_id%03d.dat",id));
      //c1->cd(1);
      tree->Draw("nebulaTURaw>>ht1(4500,0.5,4500.5)",Form("nebulaID==%d",id),"goff");
      ht1 = (TH1D*)gDirectory->Get("ht1");
      ht1->SetNameTitle(Form("NEBULA TDC_up ID%d",id),Form("NEBULA TDC_up ID%d",id));
      ht1->GetXaxis()->SetTitle("channel (ch)");
      ts1->Search(ht1,2,"goff");
      cout<<"ID"<<id<<"(u): "<<ts1->GetNPeaks()<<endl;
      Double_t *peakX1 = new Double_t[ts1->GetNPeaks()];
      Int_t *index1 = new Int_t[ts1->GetNPeaks()];
      Double_t *ch1 = new Double_t[ts1->GetNPeaks()];
      Double_t *time1 = new Double_t[ts1->GetNPeaks()];
      Double_t *res1 = new Double_t[ts1->GetNPeaks()];
      TMath::Sort(ts1->GetNPeaks(),ts1->GetPositionX(),index1,false);
      for(Int_t i=0;i<ts1->GetNPeaks();i++)
	{
	  peakX1[i] = (ts1->GetPositionX())[index1[i]];
	  time1[i] = 10*(i+1);
	  Double_t min1 = peakX1[i]-10;
	  Double_t max1 = peakX1[i]+10;
	  Double_t center1 = peakX1[i];
	  Double_t range1 = max1-min1;
	  tree->Draw(Form("nebulaTURaw>>hi%d(%lf,%lf,%lf)",i,range1,min1,max1),Form("nebulaID==%d&&abs(nebulaTURaw-%lf)<%lf",id,center1,range1),"goff");
	  hi1 = (TH1D*)gDirectory->Get(Form("hi%d",i));
	  peakX1[i] = hi1->GetMean();
	  fout1<<peakX1[i]<<" "<<time1[i]<<endl;
	  cout<<(double)i/(ts1->GetNPeaks())*100<<"% done" <<"\r";
	  cout.flush();
	}
      
      //c1->cd(5);
      tree->Draw("nebulaTDRaw>>ht2(4500,0.5,4500.5)",Form("nebulaID==%d",id),"goff");
      ht2 = (TH1D*)gDirectory->Get("ht2");
      ht2->SetNameTitle(Form("NEBULA TDC_down ID%d",id),Form("NEBULA TDC_down ID%d",id));
      ht2->GetXaxis()->SetTitle("channel (ch)");
      ts2->Search(ht2,2,"goff");
      cout<<"ID"<<id<<"(d): "<<ts2->GetNPeaks()<<endl;
      Double_t *peakX2 = new Double_t[ts2->GetNPeaks()];
      Int_t *index2 = new Int_t[ts2->GetNPeaks()];
      Double_t *ch2 = new Double_t[ts2->GetNPeaks()];
      Double_t *time2 = new Double_t[ts2->GetNPeaks()];
      Double_t *res2 = new Double_t[ts2->GetNPeaks()];
      TMath::Sort(ts2->GetNPeaks(),ts2->GetPositionX(),index2,false);
      for(Int_t i=0;i<ts2->GetNPeaks();i++)
	{
	  peakX2[i] = (ts2->GetPositionX())[index2[i]];
	  time2[i] = 10*(i+1);
	  Double_t min2 = peakX2[i]-10;
	  Double_t max2 = peakX2[i]+10;
	  Double_t center2 = peakX2[i];
	  Double_t range2 = max2-min2;
	  tree->Draw(Form("nebulaTDRaw>>hi%d(%lf,%lf,%lf)",i,range2,min2,max2),Form("nebulaID==%d&&abs(nebulaTDRaw-%lf)<%lf",id,center2,range2),"goff");
	  hi1 = (TH1D*)gDirectory->Get(Form("hi%d",i));
	  peakX2[i] = hi1->GetMean();
	  fout2<<peakX2[i]<<" "<<time2[i]<<endl;
	  cout<<(double)i/(ts2->GetNPeaks())*100<<"% done" <<"\r";
	  cout.flush();
	}

      fout1.close();
      fout2.close();
      
    }


  
}



void tdc_fit(int order=1)
{


  TF1 *fit1;TF1 *fit2;TF1 *fit3;TF1 *fit4;TF1 *fit5;TF1 *fit6;
  TF1 *func_fit;

  TH1 *hres1 = new TH1D();
  TH1 *hres2 = new TH1D();
  
  TString numberingList[6] = {"1st","2nd","3rd","4th","5th","6th"};
  TString numbering = numberingList[order-1];

  ofstream fout1(Form("./dat/NEBULA/par_%s_time_vs_channel_u.dat",numbering.Data()));
  ofstream fout2(Form("./dat/NEBULA/par_%s_time_vs_channel_d.dat",numbering.Data()));
  
  TCanvas *c2 = new TCanvas("c2","c2",1200,600);
  c2->Divide(2,1);
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  c1->Divide(3,2);

  c1->Print(Form("./fig/tcal_NEBULA_%s.pdf[",numbering.Data()));
  c2->Print(Form("./fig/tcal_NEBULA_res_%s.pdf[",numbering.Data()));
  

  
  for(Int_t id=1;id<=144;id++)
    {
      hres1->Delete();
      hres2->Delete();

      Double_t peakX1[50]={0,};
      Double_t time1[50] ={0,};
      Double_t res1[50]  ={0,};
      
      Double_t peakX2[50]={0,};
      Double_t time2[50] ={0,};
      Double_t res2[50]  ={0,};
      
      Int_t lineNum1 = 0;
      ifstream fin1(Form("./dat/NEBULA/time_vs_channel_u_id%03d.dat",id));
      while(1)
	{
	  if(fin1.fail())
	    {
	      lineNum1--;
	      break;
	    }
	  fin1>>peakX1[lineNum1]>>time1[lineNum1];
	  lineNum1++;
	}

      fin1.close();

      c1->cd(1);
      TGraph *gr1 = new TGraph(lineNum1,peakX1,time1);
      gr1->SetMarkerStyle(24);
      gr1->SetTitle(Form("NEBULA TDC_up ID%d;channel (ch);time (nsec)",id));
      gr1->Draw("AP");

      Double_t lrange = peakX1[2]-50;
      Double_t rrange = peakX1[lineNum1-5]+50;

      fit1 = new TF1("fit1","pol1",lrange,rrange);
      gr1->Fit("fit1","rq");

      
      if(order>=2)
	{
	  fit2 = new TF1("fit2","pol2",lrange,rrange);
	  fit2->SetParameter(0,fit1->GetParameter(0));
	  fit2->SetParameter(1,fit1->GetParameter(1));
	  fit2->SetParameter(2,0);
	  gr1->Fit("fit2","rq");
	}

      if(order>=3)
	{
	  fit3 = new TF1("fit3","pol3",lrange,rrange);
	  fit3->SetParameter(0,fit2->GetParameter(0));
	  fit3->SetParameter(1,fit2->GetParameter(1));
	  fit3->SetParameter(2,fit2->GetParameter(2));
	  fit3->SetParameter(3,0);
	  gr1->Fit("fit3","rq");
	}

      if(order>=4)
	{
	  fit4 = new TF1("fit4","pol4",lrange,rrange);
	  fit4->SetParameter(0,fit3->GetParameter(0));
	  fit4->SetParameter(1,fit3->GetParameter(1));
	  fit4->SetParameter(2,fit3->GetParameter(2));
	  fit4->SetParameter(3,fit3->GetParameter(3));
	  fit4->SetParameter(4,0);
	  gr1->Fit("fit4","rq");
	}

      if(order>=5)
	{
	  fit5 = new TF1("fit5","pol5",lrange,rrange);
	  fit5->SetParameter(0,fit4->GetParameter(0));
	  fit5->SetParameter(1,fit4->GetParameter(1));
	  fit5->SetParameter(2,fit4->GetParameter(2));
	  fit5->SetParameter(3,fit4->GetParameter(3));
	  fit5->SetParameter(4,fit4->GetParameter(4));
	  fit5->SetParameter(5,0);
	  gr1->Fit("fit5","rq");
	}

      if(order>=6)
	{
	  fit6 = new TF1("fit6","pol6",lrange,rrange);
	  fit6->SetParameter(0,fit5->GetParameter(0));
	  fit6->SetParameter(1,fit5->GetParameter(1));
	  fit6->SetParameter(2,fit5->GetParameter(2));
	  fit6->SetParameter(3,fit5->GetParameter(3));
	  fit6->SetParameter(4,fit5->GetParameter(4));
	  fit6->SetParameter(5,fit5->GetParameter(5));
	  fit6->SetParameter(6,0);
	  gr1->Fit("fit6","rq");
	}
      
      if(order==1)
	func_fit = fit1;
      else if(order==2)
	func_fit = fit2;
      else if(order==3)
	func_fit = fit3;
      else if(order==4)
	func_fit = fit4;
      else if(order==5)
	func_fit = fit5;
      else if(order==6)
	func_fit = fit6;

      for(int j=0;j<=order;j++)
	{
	  fout1<<func_fit->GetParameter(j)<<" ";
	}
      fout1<<endl;
      
      c1->cd(2);
      hres1 = new TH1D("hres1","hres1",15,-0.1,0.1);
      for(Int_t i=2;i<lineNum1-5;i++)
	{
	  //cout<<time1[i]<<" "<<func_fit->Eval(peakX1[i])<<endl;
	  res1[i] = time1[i] - func_fit->Eval(peakX1[i]);
	  hres1->Fill(res1[i]);
	}
      hres1->SetTitle(Form("NEBULA time residual (up) ID%d;t_{res} (nsec)",id));
      hres1->Draw();
      
      c1->cd(3);
      TGraph *gr11 = new TGraph(lineNum1,peakX1,res1);
      gr11->SetTitle("NEBULA time residual (up);channel (ch);t_{res} (nsec)");
      gr11->GetYaxis()->SetRangeUser(-0.1,0.1);
      gr11->SetMarkerStyle(24);
      gr11->Draw("AP");



      
      Int_t lineNum2 = 0;
      ifstream fin2(Form("./dat/NEBULA/time_vs_channel_d_id%03d.dat",id));
      while(1)
	{
	  if(fin2.fail())
	    {
	      lineNum2--;
	      break;
	    }
	  fin2>>peakX2[lineNum2]>>time2[lineNum2];
	  lineNum2++;
	}
      
      fin2.close();
      
      c1->cd(4);
      TGraph *gr2 = new TGraph(lineNum2,peakX2,time2);
      gr2->SetMarkerStyle(24);
      gr2->SetTitle(Form("NEBULA TDC_down ID%d;channel (ch);time (nsec)",id));
      gr2->Draw("AP");
      Double_t lrange2 = peakX2[2]-50;
      Double_t rrange2 = peakX2[lineNum2-5]+50;
      fit1 = new TF1("fit1","pol1",lrange2,rrange2);
      gr2->Fit("fit1","rq");


      if(order>=2)
	{
	  fit2 = new TF1("fit2","pol2",lrange2,rrange2);
	  fit2->SetParameter(0,fit1->GetParameter(0));
	  fit2->SetParameter(1,fit1->GetParameter(1));
	  fit2->SetParameter(2,0);
	  gr2->Fit("fit2","rq");
	}

      if(order>=3)
	{
	  fit3 = new TF1("fit3","pol3",lrange2,rrange2);
	  fit3->SetParameter(0,fit2->GetParameter(0));
	  fit3->SetParameter(1,fit2->GetParameter(1));
	  fit3->SetParameter(2,fit2->GetParameter(2));
	  fit3->SetParameter(3,0);
	  gr2->Fit("fit3","rq");
	}

      if(order>=4)
	{
	  fit4 = new TF1("fit4","pol4",lrange2,rrange2);
	  fit4->SetParameter(0,fit3->GetParameter(0));
	  fit4->SetParameter(1,fit3->GetParameter(1));
	  fit4->SetParameter(2,fit3->GetParameter(2));
	  fit4->SetParameter(3,fit3->GetParameter(3));
	  fit4->SetParameter(4,0);
	  gr2->Fit("fit4","rq");
	}

      if(order>=5)
	{
	  fit5 = new TF1("fit5","pol5",lrange2,rrange2);
	  fit5->SetParameter(0,fit4->GetParameter(0));
	  fit5->SetParameter(1,fit4->GetParameter(1));
	  fit5->SetParameter(2,fit4->GetParameter(2));
	  fit5->SetParameter(3,fit4->GetParameter(3));
	  fit5->SetParameter(4,fit4->GetParameter(4));
	  fit5->SetParameter(5,0);
	  gr2->Fit("fit5","rq");
	}

      if(order>=6)
	{
	  fit6 = new TF1("fit6","pol6",lrange2,rrange2);
	  fit6->SetParameter(0,fit5->GetParameter(0));
	  fit6->SetParameter(1,fit5->GetParameter(1));
	  fit6->SetParameter(2,fit5->GetParameter(2));
	  fit6->SetParameter(3,fit5->GetParameter(3));
	  fit6->SetParameter(4,fit5->GetParameter(4));
	  fit6->SetParameter(5,fit5->GetParameter(5));
	  fit6->SetParameter(6,0);
	  gr2->Fit("fit6","rq");
	}
      
      if(order==1)
	func_fit = fit1;
      else if(order==2)
	func_fit = fit2;
      else if(order==3)
	func_fit = fit3;
      else if(order==4)
	func_fit = fit4;
      else if(order==5)
	func_fit = fit5;
      else if(order==6)
	func_fit = fit6;


      for(int j=0;j<=order;j++)
	{
	  fout2<<func_fit->GetParameter(j)<<" ";
	}
      fout2<<endl;
      
      


      c1->cd(5);
      hres2 = new TH1D("hres2","hres2",15,-0.1,0.1);

      for(Int_t i=2;i<lineNum2-5;i++)
	{
	  //cout<<time2[i]<<" "<<func_fit->Eval(peakX2[i])<<endl;
	  res2[i] = time2[i] - func_fit->Eval(peakX2[i]);
	  hres2->Fill(res2[i]);
	}
      hres2->SetTitle(Form("NEBULA time residual (down) ID%d;t_{res} (nsec)",id));
      hres2->Draw();
      
      c1->cd(6);
      TGraph *gr22 = new TGraph(lineNum2,peakX2,res2);
      gr22->SetTitle("NEBULA time residual (down);channel (ch);t_{res} (nsec)");
      gr22->GetYaxis()->SetRangeUser(-0.1,0.1);
      gr22->SetMarkerStyle(24);
      gr22->Draw("AP");




      

      c2->cd(1);
      if(id==1)
	gr11->Draw("AL");
      else
	gr11->Draw("SAMEL");


      c2->cd(2);
      if(id==1)
	gr22->Draw("AL");
      else
	gr22->Draw("SAMEL");

      
      c2->Update();
      c1->Update();


      c1->Print(Form("./fig/tcal_NEBULA_%s.pdf",numbering.Data()));

      
      //getchar();    
      
    }

  c1->Print(Form("./fig/tcal_NEBULA_%s.pdf]",numbering.Data()));
  c2->Print(Form("./fig/tcal_NEBULA_res_%s.pdf",numbering.Data()));
  c2->Print(Form("./fig/tcal_NEBULA_res_%s.pdf]",numbering.Data()));
  

  /*
    delete fit1;
    delete fit2;
    delete fit3;
    delete fit4;
    delete fit5;
    delete fit6;
    delete func_fit;
  */
  
  fout1.close();
  fout2.close();
}

