{
  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  TCanvas *c1 = new TCanvas("c1","c1",1800,600);
  c1->Divide(3,1);
  c1->cd(1);
  chain->Draw("nebulaTOF:nebulaQPed>>h1(1000,0,3000,200,-270,-200)","nebulaID==15","colz");
  c1->cd(2);
  chain->Draw("nebulaTOFU:nebulaQPed>>h2(1000,0,3000,200,0,100)","nebulaID==15","colz");
  c1->cd(3);
  chain->Draw("nebulaTOFU:nebulaQPed>>h3(1000,0,3000,400,-100,300)","abs(nebulaID-15.5)<15&&abs(nebulaY)<10","colz");
}
