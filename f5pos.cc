{
  TFile *file = new TFile("./root/run0261.root.PPAC","r");
  TTree *tree;
  gDirectory->GetObject("PPAC",tree);
  tree->BuildIndex("RunNum","EventNum");
  tree->AddFriend("BeamPla","./root/run0261.root.BeamPla");
  tree->GetFriend("BeamPla")->BuildIndex("RunNum","EventNum");

  TFile *file2 = new TFile("./cut/f5ppac.root","r");
  TCutG *f5cut;
  gDirectory->GetObject("f5cut",f5cut);
  
  tree->SetAlias("tofc","TOF3_13-3*dT5");
  tree->SetAlias("elossc","icbEloss-0.1*tofc-40");

  TCut cut_30ne = "abs(elossc-18)<2 && abs(tofc+389.5)<1.5";

  ofstream fout("./dat/f5ppac_par.dat");
  
  TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
  c2->Divide(2,2);
  c2->cd(1);
  tree->Draw("ppac1AX:dT5Raw>>h1(400,-500,40,400,-150,150)","ppac1AX>-999&&f5cut","col");
  tree->Draw("ppac1AX:dT5Raw>>hp1(200,-500,40)","ppac1AX>-999&&f5cut","prof");
  TF1 *fit1 = new TF1("fit1","pol5",-500,40);
  hp1->Fit(fit1,"rq","",-440,10);
  h1->Draw("COLZ");
  hp1->Draw("SAME");
  fit1->Draw("SAME");

  
  c2->cd(2);
  f5cut->SetVarY("ppac1BX");
  tree->Draw("ppac1BX:dT5Raw>>h2(400,-500,40,400,-150,150)","ppac1BX>-999&&f5cut","col");
  tree->Draw("ppac1BX:dT5Raw>>hp2(400,-500,40,400,-150,150)","ppac1BX>-999&&f5cut","prof");
  TF1 *fit2 = new TF1("fit2","pol5",-500,40);
  hp2->Fit(fit2,"rq","",-440,10);
  h2->Draw("colz");
  hp2->Draw("same");
  fit2->Draw("SAMe");
  
  c2->cd(3);
  f5cut->SetVarY("ppac2AX");
  tree->Draw("ppac2AX:dT5Raw>>h3(400,-500,40,400,-150,150)","ppac2AX>-999&&f5cut","col");
  tree->Draw("ppac2AX:dT5Raw>>hp3(400,-500,40,400,-150,150)","ppac2AX>-999&&f5cut","prof");
  TF1 *fit3 = new TF1("fit3","pol5",-500,40);
  hp3->Fit(fit3,"rq","",-440,10);
  h3->Draw("colz");
  hp3->Draw("same");
  fit3->Draw("SAME");
  
  c2->cd(4);
  f5cut->SetVarY("ppac2BX");
  tree->Draw("ppac2BX:dT5Raw>>h4(400,-500,40,400,-150,150)","ppac2BX>-999&&f5cut","col");
  tree->Draw("ppac2BX:dT5Raw>>hp4(400,-500,40,400,-150,150)","ppac2BX>-999&&f5cut","prof");
  TF1 *fit4 = new TF1("fit4","pol5",-500,40);
  hp4->Fit(fit4,"rq","",-440,10);
  h4->Draw("colz");
  hp4->Draw("same");
  fit4->Draw("same");

  for(Int_t i=0;i<6;i++)
    fout<<fit1->GetParameter(i)<<" ";
  fout<<endl;

  for(Int_t i=0;i<6;i++)
    fout<<fit2->GetParameter(i)<<" ";
  fout<<endl;

  for(Int_t i=0;i<6;i++)
    fout<<fit3->GetParameter(i)<<" ";
  fout<<endl;

  for(Int_t i=0;i<6;i++)
    fout<<fit4->GetParameter(i)<<" ";
  fout<<endl;

  fout.close();
}
