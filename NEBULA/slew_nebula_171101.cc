void slew_nebula()
{
  TChain *chain = new TChain("NEBslew","NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  
  chain->Draw("nebulaTA:nebulaQPed>>h1(4000,0,4000,200,5,100)","nebulaID>0&&nebulaID<31","colz");
  
}
