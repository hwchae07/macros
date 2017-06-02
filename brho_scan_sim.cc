//void brho_scan_sim()
{
  TString filename[34] = {"brho_scan_run_1,10T.root",  "brho_scan_run_1,95T.root",  "brho_scan_run_2,55T.root",
			  "brho_scan_run_1,2T.root" ,  "brho_scan_run_1,9T.root" ,  "brho_scan_run_2,5T.root",
			  "brho_scan_run_1,30T.root",  "brho_scan_run_2,00T.root",  "brho_scan_run_2,65T.root",
			  "brho_scan_run_1,4T.root",   "brho_scan_run_2,05T.root",  "brho_scan_run_2,6T.root",
			  "brho_scan_run_1,55T.root",  "brho_scan_run_2,10T.root",  "brho_scan_run_2,75T.root",
			  "brho_scan_run_1,5T.root",   "brho_scan_run_2,15T.root",  "brho_scan_run_2,7T.root",
			  "brho_scan_run_1,65T.root",  "brho_scan_run_2,25T.root",  "brho_scan_run_2,80T.root",
			  "brho_scan_run_1,6T.root" ,  "brho_scan_run_2,2T.root" ,  "brho_scan_run_2,85T.root",
			  "brho_scan_run_1,75T.root",  "brho_scan_run_2,35T.root",  "brho_scan_run_2,95T.root",
			  "brho_scan_run_1,7T.root" ,  "brho_scan_run_2,3T.root" ,  "brho_scan_run_2,9T.root",
			  "brho_scan_run_1,85T.root",  "brho_scan_run_2,45T.root",
			  "brho_scan_run_1,8T.root" ,  "brho_scan_run_2,4T.root"};

  TString filepath0 = "./root_simtrace/";
  TString filepath;
  Int_t count=0;
  //count the data number//
  for(Int_t i=0;i<34;i++)
    {
      filepath = filepath0+filename[i];
      TFile *file = new TFile(filepath,"r");
      TTree *tr = (TTree*)file->Get("tree");
      if( tr->GetEntries("hodid>0") )
	count++;
      tr->Delete();
      file->Delete();
    }
  //count the data number//
  Double_t *hodt = new Double_t[count];
  Double_t *hodterr = new Double_t[count];
  Double_t *hodx = new Double_t[count];
  Double_t *hodxerr = new Double_t[count];

  
  //  TFile *f1 = new TFile("./root_simtrace/brho_scan_run_2,4T.root","r");

  count = 0;
  for(Int_t i=0;i<34;i++)
    {
      filepath = filepath0+filename[i];
      cout<<filepath<<endl;
      TFile *f1 = new TFile(filepath,"r");
      TTree *tree = (TTree*)f1->Get("tree");
      
      if( tree->GetEntries("hodid>0") )
	{
	  TCanvas *c1 = new TCanvas("c1","c1",1200,600);
	  c1->Divide(2,1);
	  c1->cd(1);
	  tree->Draw("hodx>>h1","hodid>0");
	  TH1 *h1 = (TH1D*)gDirectory->Get("h1");
	  Double_t centerx = h1->GetBinCenter(h1->GetMaximumBin());
	  Double_t rangex = 10;
	  TF1 *fit1 = new TF1("fit1","gaus",centerx-rangex,centerx+rangex);
	  h1->Fit(fit1,"rq","",centerx-rangex,centerx+rangex);
	  //h1->Draw();
	  hodx[count] = fit1->GetParameter(1);
	  hodxerr[count] = fit1->GetParameter(2);

	  c1->cd(2);
	  tree->Draw("hodt>>h2","hodid>0");
	  TH1 *h2 = (TH1D*)gDirectory->Get("h2");
	  Double_t centert = h2->GetBinCenter(h2->GetMaximumBin());
	  Double_t ranget = 0.06;
	  TF1 *fit2 = new TF1("fit2","gaus",centert-ranget,centert+ranget);
	  h2->Fit(fit2,"rq","",centert-ranget,centert+ranget);
	  //h2->Draw();
	  hodt[count] = fit2->GetParameter(1);
	  hodterr[count] = fit2->GetParameter(2);

	  fit1->Delete();
	  h1->Delete();
	  fit2->Delete();
	  h2->Delete();
	  c1->Delete();
	  count++;
	}
      
      tree->Delete();
      f1->Delete();
    }
  //else
  //cout<<"nothing"<<endl;

  //for(Int_t i=0;i<count;i++)
  //cout<<hodt[i]<<" "<<hodx[i]<<endl;

  TGraphErrors *gr1 = new TGraphErrors(count,hodx,hodt,hodxerr,hodterr);
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  gr1->SetMarkerStyle(24);
  gr1->Draw("AP");
  TF1 *fithod = new TF1("fithod","pol5",-1200,1200);
  fithod->SetParameter(0,5.18e+1);
  fithod->SetParameter(1,1.35e-3);
  fithod->SetParameter(2,4.68e-7);
  fithod->SetParameter(3,3.29e-12);
  fithod->SetParameter(4,1e-17);
  fithod->SetParameter(5,1e-22);
  gr1->Fit("fithod","rq","",-1200,1200);

  cout<<"ID  "<<"\t"<<"center"<<"\t"<<"HODT"<<endl;
  for(Int_t i=1;i<=24;i++)
    {
      Double_t temp = 1250 - i*100;
      cout<<"ID"<<i<<"\t"<<temp<<"\t"<<fithod->Eval(temp)<<endl;
    }
  
}
