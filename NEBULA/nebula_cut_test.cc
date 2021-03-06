//void nebula_cut_test(int i=0)
{
  /*
    Int_t numBad1[42] = {1,2,3,4,5,27,28,29,30
    ,31,32,33,34,35,36,57,58,59,60
    ,61,62,63,64,65,66,87,88,89,90
    ,91,92,93,94,95,96,97,100,102,106,118,119,120};

    //9 + 10 + 10 + 13 = 42
    */

  Int_t numBad1[37] = {1,2,3,4,5,27,28,29,30
		       ,31,32,33,34,35,36,57,58,59,60
		       ,61,62,63,64,65,66,87,88,89,90
		       ,95,96,97,100,102,106,118,119};
  //9 + 10 + 10 + 8 = 37

  
  // no gamma-ray line //
  //Int_t numBad2[] = {61,62,63,64,90,91,92,93,94,95,96,119,120};
  Int_t numBad2[5] = {91,92,93,94,120};
  // no gamma-ray line //

  
  TChain *chain = new TChain("NEBslew","NEBslew");
  //chain->Add("./root/run0146.root.NEBslew");
  //chain->Add("./root/run014[6-7].root.NEBslew");
  chain->Add("./root/run014[6-9].root.NEBslew");
  chain->Add("./root/run015[0-4].root.NEBslew");

  Int_t i = 27;
  
  Int_t id = numBad1[i];


  Int_t layer;
  if(i<9)
    layer=1;
  else if(i<19)
    layer=2;
  else if(i<29)
    layer=3;
  else if(i<37)
    layer=4;
  else
    return;

  TFile *file = new TFile(Form("./cut/nebula_slew_add_layer%d.root",layer),"r");
  TCutG *cut = (TCutG*)file->Get(Form("cutID%d",id));

  TCanvas *c1 = new TCanvas(Form("c%d",i),"c1",1200,600);
  c1->Divide(2,1);

  TH2 *h0;
  TH2 *h1;
  TH1 *h2;
  TH2 *h3;
  TProfile *h4;
  TH1 *h5;
  
  c1->cd(1);
  chain->Draw("nebulaTOF:nebulaQPed>>h0(200,30,430,200,-270,-240)",Form("nebulaID==%d",id),"colz");
  h0 = (TH2D*)gDirectory->Get("h0");
  h0->GetZaxis()->SetRangeUser(0.1,10);
  cut->Draw("SAME"); 
  
    
  //chain->Draw(">>el",Form("nebulaID==%d&&cutID%d",id,id));
  /*
    chain->Draw(">>el",Form("cutID%d",id));
    TEventList *el = (TEventList*)gDirectory->Get("el");
    chain->SetEventList(el);
  */

  /*
    chain->Draw(">>el1",Form("nebulaID==%d",id));
    TEventList *el1 = (TEventList*)gDirectory->Get("el1");
    chain->SetEventList(el1);
  
    chain->Draw(">>el2",Form("cutID%d",id));
    TEventList *el2 = (TEventList*)gDirectory->Get("el2");
    chain->SetEventList(el2);
  */
  
  c1->cd(2);
  //chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-270,-240)",Form("nebulaID==%d",id),"colz");
  //chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-270,-240)","","colz");
  chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-270,-240)",Form("cutID%d&&nebulaID==%d",id,id),"colz");
  //chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-270,-240)",Form("nebulaID==%d&&cutID%d",id,id),"colz");
  //chain->Draw("nebulaTOF:nebulaQPed>>h1(200,30,430,200,-270,-240)",Form("cutID%d",id),"colz");
  h1 = (TH2D*)gDirectory->Get("h1");
  h1->GetZaxis()->SetRangeUser(0.1,10);

}
