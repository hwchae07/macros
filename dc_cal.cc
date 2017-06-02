#include "dc_cal.hh"

void LoadTDCDistribution(const char *filename);
  
void analysis(Int_t number = 1)
{
  fBDCParameters = TArtSAMURAIParameters::Instance();
  fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC1.xml");
  fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC2.xml");

  fCalibBDC1Hit = new TArtCalibBDC1Hit;
  fCalibBDC2Hit = new TArtCalibBDC2Hit;
  
  fCalibBDC1Track = new TArtCalibBDC1Track;
  fCalibBDC2Track = new TArtCalibBDC2Track;

  LoadTDCDistribution((char*)"./db/dc/s027_run0167.root");

  
  TArtEventStore *eventStore = new TArtEventStore;
  eventStore->Open("./sdaq02/run0167.ridf");
  Int_t nEvent = 0;
  Double_t tdc[5000];
  Double_t DL[5000];

  TFile *file = new TFile("./root/tdc_test.root","recreate");
  TTree *bdc = new TTree("bdc","bdc");  
  //Int_t bdc1idNum;
  //Int_t bdc1id;
  //bdc->Branch("bdc1idNum",&bdc1idNum,"bdc1idNum/I");
  //bdc->Branch("bdc1layer",&bdc1layer,"bdc1layer[100]/I");
  //bdc->Branch("bdc1tdc",&bdc1tdc,"bdc1tdc[100]/I");
  //bdc->Branch("bdc1dl",&bdc1dl,"bdc1dl[100]/D");

  //Branch for std::vector//
  bdc->Branch("bdc1idNum",&bdc1idNum,"bdc1idNum/I");
  bdc->Branch("bdc1layer",&bdc1layer);
  bdc->Branch("bdc1tdc",&bdc1tdc);
  bdc->Branch("bdc1dl",&bdc1dl);
  //Branch for std::vector//

  
  Int_t nMaxEvent = number;
  
  while(eventStore->GetNextEvent()&&(!nMaxEvent || nEvent<nMaxEvent))
    {

      if(!(nEvent % 1000))
	{
	  cout<<"Event Number : "<<nEvent<<"\r";
	  cout.flush();
	}
      
      fCalibBDC1Hit->ClearData();
      fCalibBDC1Track->ClearData();
      fCalibBDC1Hit->ReconstructData();
      fCalibBDC1Track->ReconstructData();
      
      //fCalibBDC2Track->ClearData();
      //fCalibBDC2Track->ReconstructData();

      
      
      Int_t NumDCHit = fCalibBDC1Hit->GetNumDCHit();
      bdc1idNum = NumDCHit;

      /*
	bdc1layer = new Int_t[bdc1idNum];
	bdc1tdc = new Int_t[bdc1idNum];
	bdc1dl = new Double_t[bdc1idNum];
      */
      for(Int_t i=0;i<NumDCHit;i++)
	{
	  TArtDCHit *hit = fCalibBDC1Hit->GetDCHit(i);
	  
	  bdc1layer.push_back( hit->GetLayer() );
	  bdc1tdc.push_back( hit->GetTDC() );

	  //bdc1id = hit->GetHitID();
	  //bdc1layer[bdc1id] = hit->GetLayer();
	  //bdc1tdc[bdc1id] = hit->GetTDC();
	}
      
      Int_t NumDCTrack = fCalibBDC1Track->GetNumDCTrack();

      for(Int_t track=0;track<NumDCTrack;track++)
	{
	  TArtDCTrack *trk1 = fCalibBDC1Track->GetDCTrack(track);

	  Int_t NumHitLayer = trk1->GetNumHitLayer();
	  for(Int_t layer=0;layer<NumHitLayer;layer++)
	    {
	      for(Int_t j=0;j<NumDCHit;j++)
		{
		  if(trk1->GetHitID(layer) == j)
		    bdc1dl.push_back( trk1->GetDriftLength(layer) );
		  //bdc1dl[j] = trk1->GetDriftLength(layer);
		}
	    }
	}

      /*
	cout<<"-------------------"<<endl;
	for(Int_t j=0;j<NumDCHit;j++)
	cout<<j<<"\t"<<bdc1layer[j]<<"\t"<<bdc1tdc[j]<<"\t"<<bdc1dl[j]<<endl;
	cout<<"-------------------"<<endl;
      */

      bdc->Fill();      
      nEvent++;
           
      //delete bdc1dl;
      //delete bdc1tdc;
      //delete bdc1layer;
      bdc1layer.clear();
      bdc1tdc.clear();
      bdc1dl.clear();
      
    }

  file->Write();
  file->Close();
  
  delete eventStore;
}
  void getDL(Int_t number = 1)
  {
    fBDCParameters = TArtSAMURAIParameters::Instance();
    fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC1.xml");
    fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC2.xml");

    fCalibBDC1Hit = new TArtCalibBDC1Hit;
    fCalibBDC2Hit = new TArtCalibBDC2Hit;
  
    fCalibBDC1Track = new TArtCalibBDC1Track;
    fCalibBDC2Track = new TArtCalibBDC2Track;

    LoadTDCDistribution((char*)"./db/dc/s027_run0167.root");

  
    TArtEventStore *eventStore = new TArtEventStore;
    eventStore->Open("./sdaq02/run0167.ridf");
    Int_t nEvent = 0;
    Double_t tdc[5000];
    Double_t DL[5000];

    Int_t nMaxEvent = number;
  
    while(eventStore->GetNextEvent()&&(!nMaxEvent || nEvent<nMaxEvent))
      {
	fCalibBDC1Hit->ClearData();
	fCalibBDC1Track->ClearData();
	fCalibBDC1Hit->ReconstructData();
	fCalibBDC1Track->ReconstructData();
      
	//fCalibBDC2Track->ClearData();
	//fCalibBDC2Track->ReconstructData();


	cout<<"--------------------------------------"<<endl;
	cout<<"Event Number : "<<nEvent<<endl;

	Int_t NumDCHit = fCalibBDC1Hit->GetNumDCHit();
	cout<<"Number of Hit : "<<NumDCHit<<endl;

	for(Int_t i=0;i<NumDCHit;i++)
	  {
	    TArtDCHit *hit = fCalibBDC1Hit->GetDCHit(i);
	    cout<<"layer : "<<hit->GetLayer()<<"\tID : "<<hit->GetHitID()<<"\tTDC : "<<hit->GetTDC()<<endl;
	  }
      
	Int_t NumDCTrack = fCalibBDC1Track->GetNumDCTrack();
	cout<<"Number of Track : "<<NumDCTrack<<endl;


	for(Int_t track=0;track<NumDCTrack;track++)
	  {
	    TArtDCTrack *trk1 = fCalibBDC1Track->GetDCTrack(track);
	    Int_t NumHitLayer = trk1->GetNumHitLayer();
	    cout<<"Track Number : "<<track<<endl;
	    cout<<"Number of Hit Layer : "<<NumHitLayer<<endl;
	    for(Int_t layer=0;layer<NumHitLayer;layer++)
	      {
		cout<<"DL of layer"<<layer<<" : "<<setw(9)<<trk1->GetDriftLength(layer)<<"\t";
		cout<<"ID : "<<trk1->GetHitID(layer)<<endl;
		//cout<<"Plane ID : "<<trk1->GetHitPlaneID(layer)<<endl;
	      }
	  }
	cout<<"--------------------------------------"<<endl;

	nEvent++;
      }

    delete eventStore;
  }


  void dc_cal()
  {
    TFile *file = new TFile("./root/BDC_tdc.root");
    TTree *tree = (TTree*)file->Get("BDC");

    TCanvas *c1 = new TCanvas("c1","c1",600,600);


    tree->Draw("bdc1tdc[0]>>h1(200,1600,1800)");
	  
    
    /*
      for(Int_t id=0;id<=7;id++)
      {
      tree->Draw("bdc1tdc>>h1(200,1600,1800)",Form("bdc1layer==%d",id));
      c1->Update();
      getchar();
      }
    */
  }


  void LoadTDCDistribution(const char *filename)
  {
    char myname[128];
    TFile *fdcin = new TFile(filename,"READ");
    gROOT->cd();
    TH1F *hist = NULL;

    for(Int_t i=0;i<4;i++)
      {
	sprintf(myname,"bdc1_ftdc_corr_%d",i);
	hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionX();
	fCalibBDC1Track->SetTDCDistribution(hist,i*2);
	delete hist;
	hist = NULL;

	hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionY();
	fCalibBDC1Track->SetTDCDistribution(hist,i*2+1);
	delete hist;
	hist = NULL;

	sprintf(myname,"bdc2_ftdc_corr_%d",i);
	hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionX();
	fCalibBDC2Track->SetTDCDistribution(hist,i*2);
	delete hist;
	hist = NULL;

	hist = (TH1F*)((TH2F*)fdcin->Get(myname))->ProjectionY();
	fCalibBDC2Track->SetTDCDistribution(hist,i*2+1);
	delete hist;
	hist = NULL;
      }
  }
