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

  bdc->Branch("bdc1isX",&bdc1isX);
  bdc->Branch("bdc1wireZ",&bdc1wireZ);
  bdc->Branch("bdc1trackX",&bdc1trackX);
  bdc->Branch("bdc1trackA",&bdc1trackA);
  bdc->Branch("bdc1trackY",&bdc1trackY);
  bdc->Branch("bdc1trackB",&bdc1trackB);
  bdc->Branch("bdc1hitX",&bdc1hitX);
  bdc->Branch("bdc1hitY",&bdc1hitY);
  bdc->Branch("bdc1hitZ",&bdc1hitZ);
  bdc->Branch("bdc1calX",&bdc1calX);
  bdc->Branch("bdc1calY",&bdc1calY);
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
	  bdc1wireZ.push_back( hit->GetWireZPosition() );

	  //bdc1id = hit->GetHitID();
	  //bdc1layer[bdc1id] = hit->GetLayer();
	  //bdc1tdc[bdc1id] = hit->GetTDC();
	}
      
      Int_t NumDCTrack = fCalibBDC1Track->GetNumDCTrack();

      Int_t *numberCheck = new Int_t[NumDCHit];
      for(Int_t j=0;j<NumDCHit;j++)
	numberCheck[j] = 0;
      Int_t NumHitLayer;

      for(Int_t j=0;j<NumDCHit;j++)
	{
	  for(Int_t track=0;track<NumDCTrack;track++)
	    {
	      if(numberCheck[j])
		break;
	      TArtDCTrack *trk1 = fCalibBDC1Track->GetDCTrack(track);

	      NumHitLayer = trk1->GetNumHitLayer();
	      for(Int_t layer=0;layer<NumHitLayer;layer++)
		{

		  if(trk1->GetHitID(layer) == j)
		    {
		      //cout<<"NumDCHit"<<j<<" NumTrack"<<track<<" NumLayer"<<layer<<endl;
		      //cout<<trk1->GetHitID(layer)<<" "<<j<<endl;
		      numberCheck[j] = 1;		      
		      bdc1hitZ.push_back( trk1->GetHitZPosition(layer) );
		      bdc1trackX.push_back( trk1->GetPosition(0) );
		      bdc1trackY.push_back( trk1->GetPosition(1) );
		      bdc1trackA.push_back( trk1->GetAngle(0) );
		      bdc1trackB.push_back( trk1->GetAngle(1) );
		      bdc1dl.push_back( trk1->GetDriftLength(layer) );
		      if(trk1->GetPosition(0) > -9999)
			{
			  bdc1isX.push_back(1);
			  bdc1hitX.push_back( trk1->GetHitXPosition(layer) );
			  bdc1hitY.push_back( -9999. );
			  bdc1calX.push_back( trk1->GetPosition(0)+ trk1->GetHitZPosition(layer)*trk1->GetAngle(0) );
			  bdc1calY.push_back( -9999. );
			}
		      else
			{
			  bdc1isX.push_back(0);
			  bdc1hitX.push_back( -9999. );
			  bdc1hitY.push_back( trk1->GetHitXPosition(layer) );
			  bdc1calX.push_back( -9999. );
			  bdc1calY.push_back( trk1->GetPosition(1)+ trk1->GetHitZPosition(layer)*trk1->GetAngle(1) );
			}

		      if(numberCheck[j])
			break;
		      
		    }
		
		}
	    }
	}

      /*
	cout<<"----------------------------------------------"<<endl;
	cout<<"Event Number : "<<nEvent<<endl;
	cout<<"ID\tLayer\tTDC\tDL\t\t";
	cout<<"is X?\t";
	cout<<"hitX\t\t";
	cout<<"calX\t\t";
	cout<<"hitY\t\t";
	cout<<"calY\t\t";
	cout<<endl;
	for(Int_t j=0;j<NumDCHit;j++)
	{
	cout<<j<<"\t"<<bdc1layer[j]<<"\t"<<bdc1tdc[j]<<"\t"<<setw(10)<<bdc1dl[j];
	cout<<"\t"<<bdc1isX[j];
	cout<<"\t"<<setw(8)<<bdc1hitX[j];
	cout<<"\t"<<setw(8)<<bdc1calX[j];
	cout<<"\t"<<setw(8)<<bdc1hitY[j];
	cout<<"\t"<<setw(8)<<bdc1calY[j];
	cout<<endl;
	}
	cout<<"----------------------------------------------"<<endl;
      */

      bdc->Fill();      
      nEvent++;
           
      //delete bdc1dl;
      //delete bdc1tdc;
      //delete bdc1layer;
      bdc1layer.clear();
      bdc1tdc.clear();
      bdc1dl.clear();

      bdc1trackX.clear();
      bdc1trackY.clear();
      bdc1trackA.clear();
      bdc1trackB.clear();

      bdc1hitX.clear();
      bdc1hitY.clear();
      bdc1hitZ.clear();
      
      bdc1calX.clear();
      bdc1calY.clear();
      bdc1wireZ.clear();

      bdc1isX.clear();
      
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
	  cout<<"layer : "<<hit->GetLayer()<<"\tID : "<<hit->GetHitID()<<"\tTDC : "<<hit->GetTDC();
	  cout<<"\tZ_wire : "<<hit->GetWireZPosition()<<endl;
	}
      
      Int_t NumDCTrack = fCalibBDC1Track->GetNumDCTrack();
      cout<<"Number of Track : "<<NumDCTrack<<endl;


      for(Int_t track=0;track<NumDCTrack;track++)
	{
	  TArtDCTrack *trk1 = fCalibBDC1Track->GetDCTrack(track);
	  Int_t NumHitLayer = trk1->GetNumHitLayer();

	  
	  cout<<"Track Number : "<<track<<endl;
	  cout<<"Number of Hit Layer : "<<NumHitLayer<<endl;
	  cout<<"ID\t";
	  //cout<<"layer\t";
	  cout<<"DL\t\t";
	  cout<<"trackX(Y)\t";
	  //cout<<"trackY\t\t";
	  //cout<<"trackA\t\t";
	  cout<<"hitX(Y)\t\t";
	  //cout<<"hitZ\t";
	  cout<<"calX(Y)\t\t";
	  cout<<"hit-track\t";
	  cout<<"cal-track\t";
	  //cout<<"calY\t\t";
	  cout<<endl;
	  for(Int_t layer=0;layer<NumHitLayer;layer++)
	    {
	      cout<<trk1->GetHitID(layer);
	      //cout<<"\t"<<layer;
	      cout<<"\t"<<setw(8)<<trk1->GetDriftLength(layer);
	      if(trk1->GetPosition(0)>-9000)
		{
		  cout<<"\t"<<setw(8)<<trk1->GetPosition(0);
		  cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer);
		  cout<<"\t"<<setw(8)<<trk1->GetPosition(0)+trk1->GetAngle(0)*trk1->GetHitZPosition(layer);
		  cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer)-trk1->GetPosition(0);
		  cout<<"\t"<<setw(8)<<trk1->GetAngle(0)*trk1->GetHitZPosition(layer);
		}
	      else
		{
		  cout<<"\t"<<setw(8)<<trk1->GetPosition(1);
		  cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer);
		  cout<<"\t"<<setw(8)<<trk1->GetPosition(1)+trk1->GetAngle(1)*trk1->GetHitZPosition(layer);
		  cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer)-trk1->GetPosition(1);
		  cout<<"\t"<<setw(8)<<trk1->GetAngle(1)*trk1->GetHitZPosition(layer);
		}
	      
	      //cout<<"\t"<<setw(8)<<trk1->GetAngle(0);
	      //cout<<"\t"<<setw(5)<<trk1->GetHitZPosition(layer);
	      

	      cout<<endl;
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
  TFile *file = new TFile("./root/tdc_test.root");
  TTree *tree = (TTree*)file->Get("bdc");

  tree->Draw("bdc1calX-bdc1trackX:bdc1tdc","bdc1layer==1");
    
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
