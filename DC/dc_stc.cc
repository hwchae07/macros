{
  TFile *file = new TFile("./root/tdc_test.root");
  TTree *tree = (TTree*)file->Get("bdc");



  double cut1x[6],cut1y[6];
  
  cut1x[0]=1738.84; cut1y[0]=0.0257537;
  cut1x[1]=1727.97; cut1y[1]=0.682946;
  cut1x[2]=1697.33; cut1y[2]=2.03266;
  cut1x[3]=1666.86; cut1y[3]=1.538;
  cut1x[4]=1664.01; cut1y[4]=-0.00251259;
  cut1x[5]=1738.84; cut1y[5]=0.0257537;

  //for layer 3

  
  //(layer%4) <= 1 : X
  //(layer%4) >= 2 : Y
  //0,1 : X
  //2,3 : Y
  //4,5 : X
  //6,7 : Y
  int layer = 3;
  TString layer_xy;
  int layerIsX;

  if(layer%4 < 2)
    {
      layer_xy = "bdc1calX-bdc1wireX";
      layerIsX = 1;
    }
  else
    {
      layer_xy = "bdc1calY-bdc1wireY";
      layerIsX = 0;
    }
  
  TCutG *cut1 = new TCutG("cut1",6,cut1x,cut1y);
  cut1->SetVarX("bdc1tdc");
  cut1->SetVarY(Form("%s",layer_xy.Data()));

  
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(4,4);

  
  for(int i=0;i<=15;i++)
    {
      c1->cd(i+1);
      tree->Draw(Form("abs(%s):bdc1tdc>>h%d(120,1640,1760,200,-0.5,8)",layer_xy.Data(),i),Form("bdc1layer==%d&&bdc1isX==%d&&bdc1wireID==%d",layer,layerIsX,i),"colz");
      c1->Update();
    }

  TCanvas *c2 = new TCanvas("c2","c2",600,600);
  tree->Draw(Form("abs(%s):bdc1tdc>>h2d1(120,1640,1760,200,-0.5,8)",layer_xy.Data()),Form("bdc1layer==%d&&bdc1isX==%d",layer,layerIsX),"colz");


  /*
  TCanvas *c3 = new TCanvas("c3","c3",600,600);
  tree->Draw(Form("abs(%s):bdc1tdc>>h2d2(120,1640,1760,200,-0.5,4)",layer_xy.Data()),Form("bdc1layer==%d&&bdc1isX==%d&&cut1",layer,layerIsX),"colz");
  */
  /*
  for(int i=0;i<=15;i++)
    {
      double tot = tree->GetEntries(Form("bdc1layer==%d&&bdc1isX==%d&&bdc1wireID==%d",layer,layerIsX,i));
      double cut = tree->GetEntries(Form("bdc1layer==%d&&bdc1isX==%d&&bdc1wireID==%d&&cut1",layer,layerIsX,i));
      cout<<setw(2)<<i<<" : "<<setw(6)<<tot<<" "<<setw(3)<<cut<<" "<<cut/tot*100<<endl;

    }
  */
  
}
