{
  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  TFile *file = new TFile("./root/brho_hodx.root","r");
  TTree *tree = (TTree*)file->Get("tree");

  tree->Draw("brho:hodx>>hprof","brho>-9999&&hodx>-9999","prof");
  hprof->SetMinimum(5);
  hprof->SetTitle("B#rho VS x_{HODF};x (mm);B#rho (Tm)");
  TF1 *f1 = new TF1("f1","pol5",-1200,1200);
  hprof->Fit(f1,"r","",-1200,1200);

  Double_t hod_center[25]={0,};
  for(Int_t id=1;id<=24;id++)
    {
      hod_center[id] = 1250-100*id;


      //cout<<"ID"<<id<<"\t"<<hod_center[id]<<"\t"<<f1->Eval(hod_center[id])<<endl;
      cout<<id<<" : "<<f1->Eval(hod_center[id])<<endl;
    }
}
