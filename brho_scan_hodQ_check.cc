{
  TChain *chain1 = new TChain("HOD","HOD");

  chain1->Add("./root/run013[3-8].root.HOD");

  TCanvas *c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2,1);
  c1->cd(1);
  chain1->Draw("hodQA:hodID>>h1(24,0.5,24.5,500,0,1000)","","colz");
  c1->cd(2);
  chain1->Draw("hodQ:hodID>>h2(24,0.5,24.5,500,0,1000)","","colz");
}
