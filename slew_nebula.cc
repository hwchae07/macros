
{
  //22Ne beam
  //Al target
  //gamma ray production



  TChain *chain = new TChain("NEBULA","NEBULA");
  //chain->Add("./root/run0146.root.NEBULA");
  chain->Add("./root/run014[6-9].root.NEBULA");
  chain->Add("./root/run015[0-4].root.NEBULA");


  /*
    TFile *file = new TFile("./root/run0146.root.BeamPID","r");
    TTree *tree = (TTree*)file->Get("BeamPID");
  */
  //chain->Draw("beamZ:tofc713>>h1(400,203,207,400,2,14)","","colz");
  Double_t b22ne_x[9];
  Double_t b22ne_y[9];
  b22ne_x[0]=203.995; b22ne_y[0]=10.8049;
  b22ne_x[1]=203.748; b22ne_y[1]=10.393;
  b22ne_x[2]=203.737; b22ne_y[2]=9.87246;
  b22ne_x[3]=204.08; b22ne_y[3]=9.6036;
  b22ne_x[4]=204.47; b22ne_y[4]=9.6036;
  b22ne_x[5]=204.67; b22ne_y[5]=10.0727;
  b22ne_x[6]=204.642; b22ne_y[6]=10.7248;
  b22ne_x[7]=204.431; b22ne_y[7]=10.8106;
  b22ne_x[8]=203.995; b22ne_y[8]=10.8049;

  TCutG *b22ne = new TCutG("b22ne",9,b22ne_x,b22ne_y);
  b22ne->SetVarX("tofc713");
  b22ne->SetVarY("beamZ");

  //chain->Draw("beamZ:tofc713>>h2(400,203,207,400,2,14)","b22ne","colz");

  //TH2 *h3;
  //chain->Draw("nebulaTU:nebulaID>>h3(144,0.5,144.5,300,-50,250)","","goff");
  //h3 = (TH2D*)gDirectory->Get("h3");

  TCanvas *c1 = new TCanvas("c1","c1",1500,900);
  c1->Divide(5,3);
  c1->Print("./fig/slew_check.pdf[");
  for(Int_t id=1;id<=144;id++)
    {
      c1->cd((id-1)%15+1);
      chain->Draw(Form("nebulaTU:nebulaQUPed>>hi%d(3000,0,3000,200,0,200)",id),Form("nebulaID==%d",id),"colz");
      c1->Update();
      //getchar();
      if(!(id%15)||id==144)
	c1->Print("./fig/slew_check.pdf");
    }
  c1->Print("./fig/slew_check.pdf]");
}
