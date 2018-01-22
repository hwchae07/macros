{

  //To calculate the relative HOD timing offset
  
  //gStyle->SetOptStat(11);
  
  //beam cut 22Ne, 21F and fragment cut 22Ne//
  TFile *f_cut = new TFile("./cut/brho_scan_cut.root","r");
  TCutG *b22ne, *b21f, *f22ne;
  f_cut->GetObject("b22ne",b22ne);
  f_cut->GetObject("b21f",b21f);
  f_cut->GetObject("f22ne",f22ne);
  //beam cut 22Ne, 21F and fragment cut 22Ne//

  ofstream file1("./dat/hod_relative_offset.dat");
  
  Int_t runNum = 133;
  TH1 *h1,*h2;
  
  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  c1->Divide(2,1);
  c1->Print("./fig/hod_brho.pdf[");

  for(Int_t id = 1 ; id <= 23 ; id++)
    {
      
      TH1 *hsum1 = new TH1D("hsum1","hsum1",1000,-235,-215);
      TH1 *hsum2 = new TH1D("hsum2","hsum2",1000,-235,-215);
  
      for(runNum=133;runNum<=138;runNum++)
	{
	  TString filename = Form("run%04d.root",runNum);


	  TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
	  TTree *tree = (TTree*)file->Get("BeamPla");
	  tree->BuildIndex("RunNum","EventNum");
	  tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
	  tree->GetFriend("HOD")->BuildIndex("RunNum","EventNum");
	  tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
	  tree->GetFriend("FDC")->BuildIndex("RunNum","EventNum");

	  
	  tree->Draw(">>el",Form("hodNum==2&&hodID==%d",id));
	  TEventList *el = (TEventList*)gDirectory->Get("el");
	  tree->SetEventList(el);

	  tree->Draw(">>el2",Form("hodID==%d",id+1));
	  TEventList *el2 = (TEventList*)gDirectory->Get("el2");
	  tree->SetEventList(el2);
	  
      
      
      	  tree->SetAlias("fne","abs(hodQA-450)<30");
    

	  c1->cd(1);
	  tree->Draw("hodTA-T13>>h1(50,-235,-215)",Form("b22ne&&hodID==%d",id));
	  h1 = (TH1D*)gDirectory->Get("h1");
	  hsum1->Add(h1);
	  hsum1->SetTitle(Form("HODF ID%d;TOF (arbitrary)",id));
	  if(runNum==138)
	    {
	      hsum1->Draw();
      
	      Double_t center1 = hsum1->GetBinCenter(hsum1->GetMaximumBin());
	      Double_t rms1 = hsum1->GetRMS();
	      TF1 *fit1 = new TF1("fit1","gaus");
	      fit1->SetParameter(0,hsum1->GetMaximum());
	      fit1->SetParameter(1,center1);
	      fit1->SetParameter(2,rms1);
	      hsum1->Fit("fit1","RQ","",center1-rms1,center1+rms1);
	      

	      fit1->SetParameters(fit1->GetParameters());
	      center1 = fit1->GetParameter(1);
	      rms1 = 1.5;//fit1->GetParameter(2);
	      hsum1->Fit("fit1","RQ","",center1-rms1,center1+rms1);

	      file1<<id<<","<<id+1<<"\t"<<fit1->GetParameter(1)<<"\t";

	      c1->Update();
	    }


	  c1->cd(2);
	  tree->Draw("hodTA-T13>>h2(50,-235,-215)",Form("b22ne&&hodID==%d",id+1));
	  h2 = (TH1D*)gDirectory->Get("h2");
	  hsum2->Add(h2);
	  hsum2->SetTitle(Form("HODF ID%d;TOF (arbitrary)",id+1));

	  if(runNum==138)
	    {
	      hsum2->Draw();
	  
	      Double_t center2 = hsum2->GetBinCenter(hsum2->GetMaximumBin());
	      Double_t rms2 = hsum2->GetRMS();
	      TF1 *fit2 = new TF1("fit2","gaus");
	      fit2->SetParameter(0,hsum2->GetMaximum());
	      fit2->SetParameter(1,center2);
	      fit2->SetParameter(2,rms2);
	      hsum2->Fit("fit2","RQ","",center2-rms2,center2+rms2);

	      fit2->SetParameters(fit2->GetParameters());
	      center2 = fit2->GetParameter(1);
	      rms2 = 1.5;//fit2->GetParameter(2);
	      hsum2->Fit("fit2","RQ","",center2-rms2,center2+rms2);

	      file1<<fit2->GetParameter(1)<<endl;
	      
	      c1->Update();
	    }
	  //cout<<tree->GetEntries("hodID[0]==11&&hodID[1]==12")<<endl;
	  //cout<<tree->GetEntries("hodID[0]<hodID[1]")<<endl;  
	}

      
      
      c1->Print("./fig/hod_brho.pdf");
      hsum1->Delete();
      hsum2->Delete();
    }

  c1->Print("./fig/hod_brho.pdf]");
  file1.close();


}
