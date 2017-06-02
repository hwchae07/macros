{
  TFile *file = new TFile("./root/run0275.root.BeamPID","r");
  TTree *tree = (TTree*)file->Get("BeamPID");

  TCut cut_na = "abs(beamZ-11.1)<0.5";

  tree->Draw(">>evtlist",cut_na);
  TEventList *evtlist = (TEventList*)gDirectory->Get("evtlist");
  tree->SetEventList(evtlist);

  Double_t dispersion = 3280;

  TH1 *h2;
  TF1 *fit2;

  cout<<"(x|d)\tHeight\t\tSigma\t\tHeight/Sigma"<<endl;
  
  for(Int_t i = 0 ; i<21 ; i=i+1)
    {
      dispersion = 3200 + i*20;
      //dispersion = 3279. + i*0.1;
      cout.precision(1);
      cout<<fixed<<dispersion<<"\t";
      
      /*
	TCanvas *c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	tree->Draw(Form("brho0*(1+f5X/%lf)*299.792458/beamMomU:tofc713>>h1(400,200,220,400,2.75,3.35)",dispersion),"","colz");
	c1->cd(2);*/
      tree->Draw(Form("brho0*(1+f5X/%lf)*299.792458/beamMomU>>h2(400,2.75,3.35)",dispersion),"","goff");
      h2 = (TH1D*)gDirectory->Get("h2");
      h2->Fit("gaus","rq0","goff",3.05,3.15);
      fit2 = (TF1*)h2->GetFunction("gaus");
      h2->Fit("gaus","rq0","goff",fit2->GetParameter(1)-2*fit2->GetParameter(2),fit2->GetParameter(1)+2*fit2->GetParameter(2));
      fit2 = (TF1*)h2->GetFunction("gaus");
      h2->Fit("gaus","rq0","goff",fit2->GetParameter(1)-2*fit2->GetParameter(2),fit2->GetParameter(1)+2*fit2->GetParameter(2));
      fit2 = (TF1*)h2->GetFunction("gaus");
      

      cout.precision(4);
      cout<<scientific<<fit2->GetParameter(0)<<"\t"<<fit2->GetParameter(2)<<"\t"<<fit2->GetParameter(0)/fit2->GetParameter(2)<<endl;
      //c1->Update();
      //getchar();
      //delete c1;
      /*
	cout<<"Height : " <<fit2->GetParameter(0)<<endl;
	cout<<"Mean : " <<fit2->GetParameter(1)<<endl;
	cout<<"Sigma : " <<fit2->GetParameter(2)<<endl;
	cout<<"Height/Sigma : " <<fit2->GetParameter(0)/fit2->GetParameter(2)<<endl;
      */
    }
}
