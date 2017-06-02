
bool flagEnd;

void SigHandler(int param)
{
  flagEnd = true;
}

/*
void LoadTDCDistribution(const char *filename)
{
  char myname[128];
  TFile *fdcin = new TFile(filename,"READ");
  gROOT->cd();
  TF1F *hist = NULL;

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
*/
void analysis()
{
  flagEnd = false;
  
  TArtSAMURAIParameters *fBDCParameters;
  TArtCalibBDC1Hit *fCalibBDC1Hit;
  TArtCalibBDC2Hit *fCalibBDC2Hit;
  fBDCParameters = TArtSAMURAIParameters::Instance();
  fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC1.xml");
  fBDCParameters->LoadParameter((char*)"db/SAMURAIBDC2.xml");

  fCalibBDC1Hit = new TArtCalibBDC1Hit;
  fCalibBDC2Hit = new TArtCalibBDC2Hit;

  /*
  //2017.05.19//
  TArtCalibBDC1Track *fCalibBDC1Track;
  TArtCalibBDC2Track *fCalibBDC2Track;

  fCalibBDC1Track = new TArtCalibBDC1Track;
  fCalibBDC2Track = new TArtCalibBDC2Track;

  LoadTDCDistribution((char*)"./db/dc/s027_run0167.root");
  //2017.05.19//
  */
  

  TArtEventStore *eventStore = new TArtEventStore;
  eventStore->Open("sdaq02/run0167.ridf");
  Int_t nEvent = 0;

  TFile *file = new TFile("./root/BDC_tdc.root","recreate");
  TTree *bdc = new TTree("BDC","BDC");
  
  Int_t wireID1,layer1;
  Int_t tdc1[8];
  
  bdc->Branch("bdc1wireID",&wireID1,"bdc1wireID/I");
  bdc->Branch("bdc1layer",&layer1,"bdc1layer/I");
  bdc->Branch("bdc1tdc[8]",&tdc1,"bdc1tdc[8]/I");

  Int_t wireID2,layer2;
  Int_t tdc2[8];
  bdc->Branch("bdc2wireID",&wireID2,"bdc2wireID/I");
  bdc->Branch("bdc2layer",&layer2,"bdc2layer/I");
  bdc->Branch("bdc2tdc[8]",&tdc2,"bdc2tdc[8]/I");
  
  while(!flagEnd &&eventStore->GetNextEvent())
    {
      signal(SIGINT, &SigHandler);
      if(!(nEvent %1000))
	{
	  cout<<"Event number : " <<nEvent<<"\r";
	  cout.flush();
	}
      
      fCalibBDC1Hit->ClearData();
      fCalibBDC2Hit->ClearData();
      fCalibBDC1Hit->ReconstructData();
      fCalibBDC2Hit->ReconstructData();
      /*
      fCalibBDC1Track->ClearData();
      fCalibBDC2Track->ClearData();
      fCalibBDC1Track->ReconstructionData();
      fCalibBDC2Track->ReconstructionData();
      */
      for(Int_t i=0;i<fCalibBDC1Hit->GetNumDCHit();i++)
	{
	  TArtDCHit *hit = fCalibBDC1Hit->GetDCHit(i);
	  if(hit)
	    {
	      wireID1 = hit->GetWireID();
	      layer1 = hit->GetLayer();
	      tdc1[layer1]= hit->GetTDC();
	    }
	  bdc->Fill();
	}


      for(Int_t i=0;i<fCalibBDC2Hit->GetNumDCHit();i++)
	{
	  TArtDCHit *hit = fCalibBDC2Hit->GetDCHit(i);
	  if(hit)
	    {
	      wireID2 = hit->GetWireID();
	      layer2 = hit->GetLayer();
	      tdc2[layer2] = hit->GetTDC();
	    }
	  bdc->Fill();
	}
      
      nEvent++;
      
    }
  file->Write();
  file->Close();

  delete eventStore;
  delete fCalibBDC2Hit;
  delete fCalibBDC1Hit;

}
