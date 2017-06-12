#include "dc_cal.hh"

void LoadTDCDistribution(const char *filename);

void sortid()
{
  TFile *file = new TFile("./root/tdc_test.root","r");
  TTree *tree = (TTree*)file->Get("bdc");
  Int_t NumDCHit;
  tree->SetBranchAddress("eventNum",&eventNum);
  tree->SetBranchAddress("bdc1idmax",&NumDCHit);
  tree->SetBranchAddress("bdc1id",&bdc1id);
  tree->SetBranchAddress("bdc1layer",&bdc1layer);
  tree->SetBranchAddress("bdc1tdc",&bdc1tdc);
  tree->SetBranchAddress("bdc1dl",&bdc1dl);
  tree->SetBranchAddress("bdc1isX",&bdc1isX);
  tree->SetBranchAddress("bdc1trackNum",&bdc1trackNum);
  tree->SetBranchAddress("bdc1chi2",&bdc1chi);
  tree->SetBranchAddress("bdc1wireZ",&bdc1wireZ);
  tree->SetBranchAddress("bdc1trackX",&bdc1trackX);
  tree->SetBranchAddress("bdc1trackA",&bdc1trackA);
  tree->SetBranchAddress("bdc1trackY",&bdc1trackY);
  tree->SetBranchAddress("bdc1trackB",&bdc1trackB);
  tree->SetBranchAddress("bdc1hitX",&bdc1hitX);
  tree->SetBranchAddress("bdc1hitY",&bdc1hitY);
  tree->SetBranchAddress("bdc1hitZ",&bdc1hitZ);
  tree->SetBranchAddress("bdc1calX",&bdc1calX);
  tree->SetBranchAddress("bdc1calY",&bdc1calY);
  
  TFile *file2 = new TFile("./root/tdc_sort.root","recreate");
  TTree *tree2 = new TTree("bdc","bdc");
  int eventNum2, bdc1id2, bdc1layer2, bdc1tdc2;
  int bdc1isX2, bdc1trackNum2;
  double bdc1dl2, bdc1chi2, bdc1wireZ2;
  double bdc1trackX2, bdc1trackA2, bdc1trackY2, bdc1trackB2;
  double bdc1hitX2, bdc1hitY2, bdc1hitZ2, bdc1calX2, bdc1calY2;
  tree2->Branch("eventNum",&eventNum2);
  tree2->Branch("bdc1id",&bdc1id2);
  tree2->Branch("bdc1layer",&bdc1layer2);
  tree2->Branch("bdc1tdc",&bdc1tdc2);
  tree2->Branch("bdc1dl",&bdc1dl2);
  tree2->Branch("bdc1isX",&bdc1isX2);
  tree2->Branch("bdc1trackNum",&bdc1trackNum2);
  tree2->Branch("bdc1chi2",&bdc1chi2);
  tree2->Branch("bdc1wireZ",&bdc1wireZ2);
  tree2->Branch("bdc1trackX",&bdc1trackX2);
  tree2->Branch("bdc1trackA",&bdc1trackA2);
  tree2->Branch("bdc1trackY",&bdc1trackY2);
  tree2->Branch("bdc1trackB",&bdc1trackB2);
  tree2->Branch("bdc1hitX",&bdc1hitX2);
  tree2->Branch("bdc1hitY",&bdc1hitY2);
  tree2->Branch("bdc1hitZ",&bdc1hitZ2);
  tree2->Branch("bdc1calX",&bdc1calX2);
  tree2->Branch("bdc1calY",&bdc1calY2);
  
  //for(int i = 0 ; i<3 ; i++)
  for(int i = 0 ; i<tree->GetEntries() ; i++)
    {
      /*
	cout<<endl;
	cout<<"Event "<<i<<endl;
      */
      tree->GetEntry(i);
      double chi2temp = 9999;
      int id=0;
      for(Int_t j=0;j<bdc1isX->size();j++)
	{
	  if(bdc1id->at(j) == id && bdc1chi->at(j)<chi2temp)
	    {
	      //cout<<j<<" "<<bdc1id->at(j)<<" "<<bdc1chi->at(j)<<endl;
	      chi2temp = bdc1chi->at(j);

	      eventNum2 = eventNum;
	      bdc1id2 = bdc1id->at(j);
	      bdc1layer2 = bdc1layer->at(j);
	      bdc1tdc2 = bdc1tdc->at(j);
	      bdc1dl2 = bdc1dl->at(j);
	      bdc1isX2 = bdc1isX->at(j);
	      bdc1trackNum2 = bdc1trackNum->at(j);
	      bdc1chi2 = bdc1chi->at(j);
	      bdc1wireZ2 = bdc1wireZ->at(j);
	      bdc1trackX2 = bdc1trackX->at(j);
	      bdc1trackA2 = bdc1trackA->at(j);
	      bdc1trackY2 = bdc1trackY->at(j);
	      bdc1trackB2 = bdc1trackB->at(j);
	      bdc1hitX2 = bdc1hitX->at(j);
	      bdc1hitY2 = bdc1hitY->at(j);
	      bdc1hitZ2 = bdc1hitZ->at(j);
	      bdc1calX2 = bdc1calX->at(j);
	      bdc1calY2 = bdc1calY->at(j);

	      if(j == bdc1isX->size()-1)
		tree2->Fill();
	    }
	  else if(bdc1id->at(j)>id)
	    {
	      tree2->Fill();
	      id++;
	      chi2temp = 9999;
	      j--;
	    }

	}

    }

  file2->Write();
  file2->Close();

}

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

    Int_t NumDCHit;
  
    //Branch for std::vector//
    bdc->Branch("eventNum",&eventNum);
    bdc->Branch("bdc1idmax",&NumDCHit);
    bdc->Branch("bdc1id",&bdc1id);
    bdc->Branch("bdc1layer",&bdc1layer);
    bdc->Branch("bdc1tdc",&bdc1tdc);
    bdc->Branch("bdc1dl",&bdc1dl);

    bdc->Branch("bdc1isX",&bdc1isX);
    bdc->Branch("bdc1trackNum",&bdc1trackNum);
    bdc->Branch("bdc1chi2",&bdc1chi);
  
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

	eventNum = nEvent;
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

      
      
	NumDCHit = fCalibBDC1Hit->GetNumDCHit();
	bdc1idmax = NumDCHit;
      
	Int_t NumDCTrack = fCalibBDC1Track->GetNumDCTrack();

	Int_t *numberCheck = new Int_t[NumDCHit];
	for(Int_t j=0;j<NumDCHit;j++)
	  numberCheck[j] = 0;
	Int_t NumHitLayer;

	for(Int_t j=0;j<NumDCHit;j++)
	  {
	    TArtDCHit *hit = fCalibBDC1Hit->GetDCHit(j);
	  
	    for(Int_t track=0;track<NumDCTrack;track++)
	      {
		TArtDCTrack *trk1 = fCalibBDC1Track->GetDCTrack(track);
	      
		NumHitLayer = trk1->GetNumHitLayer();

		for(Int_t layer=0;layer<NumHitLayer;layer++)
		  {
		  

		    if(trk1->GetHitID(layer) == j)
		      {
			//cout<<"NumDCHit"<<j<<" NumTrack"<<track<<" NumLayer"<<layer<<endl;
			//cout<<trk1->GetHitID(layer)<<" "<<j<<endl;
			//cout<<layer<<"\t"<<trk1->GetChi2()<<"\t"<<chi2temp<<endl;
			numberCheck[j] = 1;

			bdc1chi->push_back( trk1->GetChi2() );
			bdc1id->push_back( trk1->GetHitID(layer) );
		      
			bdc1hitZ->push_back( trk1->GetHitZPosition(layer) );
			bdc1trackNum->push_back( track );
			bdc1trackX->push_back( trk1->GetPosition(0) );
			bdc1trackY->push_back( trk1->GetPosition(1) );
			bdc1trackA->push_back( trk1->GetAngle(0) );
			bdc1trackB->push_back( trk1->GetAngle(1) );
			bdc1dl->push_back( trk1->GetDriftLength(layer) );
			if(trk1->GetPosition(0) > -9999)
			  {
			    bdc1isX->push_back(1);
			    bdc1hitX->push_back( trk1->GetHitXPosition(layer) );
			    bdc1hitY->push_back( -9999. );
			    bdc1calX->push_back( trk1->GetPosition(0)+ trk1->GetHitZPosition(layer)*trk1->GetAngle(0) );
			    bdc1calY->push_back( -9999. );
			  }
			else
			  {
			    bdc1isX->push_back(0);
			    bdc1hitX->push_back( -9999. );
			    bdc1hitY->push_back( trk1->GetHitXPosition(layer) );
			    bdc1calX->push_back( -9999. );
			    bdc1calY->push_back( trk1->GetPosition(1)+ trk1->GetHitZPosition(layer)*trk1->GetAngle(1) );
			  }

			bdc1layer->push_back( hit->GetLayer() );
			bdc1tdc->push_back( hit->GetTDC() );
			bdc1wireZ->push_back( hit->GetWireZPosition() );
		      
			//bdc->Fill();
		      
		      }
		
		  }
	      }
	  }

	/*      
	cout<<"----------------------------------------------"<<endl;
	cout<<"Event Number : "<<nEvent<<endl;
	cout<<"ID\tLayer\tTDC\tDL\t\t";
	cout<<"is X?\t";
	cout<<"trackX\t\t";
	cout<<"hitX\t\t";
	cout<<"calX\t\t";
	cout<<"hitY\t\t";
	cout<<"calY\t\t";
	cout<<endl;
	for(Int_t j=0;j<NumDCHit;j++)
	  {
	    cout<<j<<"\t"<<bdc1layer->at(j)<<"\t"<<bdc1tdc->at(j)<<"\t"<<setw(10)<<bdc1dl->at(j);
	    cout<<"\t"<<bdc1isX->at(j);
	    cout<<"\t"<<setw(8)<<bdc1trackX->at(j);
	    cout<<"\t"<<setw(8)<<bdc1hitX->at(j);
	    cout<<"\t"<<setw(8)<<bdc1calX->at(j);
	    cout<<"\t"<<setw(8)<<bdc1hitY->at(j);
	    cout<<"\t"<<setw(8)<<bdc1calY->at(j);
	    cout<<endl;
	  }
	cout<<"----------------------------------------------"<<endl;
	*/

	bdc->Fill();      
	nEvent++;
           
	bdc1layer->clear();
	bdc1tdc->clear();
	bdc1dl->clear();

	bdc1trackX->clear();
	bdc1trackY->clear();
	bdc1trackA->clear();
	bdc1trackB->clear();

	bdc1hitX->clear();
	bdc1hitY->clear();
	bdc1hitZ->clear();
      
	bdc1calX->clear();
	bdc1calY->clear();
	bdc1wireZ->clear();

	bdc1isX->clear();
	bdc1chi->clear();
	bdc1trackNum->clear();
	bdc1id->clear();
      
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
	    cout<<"\tPos_wire : "<<setw(6)<<hit->GetWirePosition();
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
	    //cout<<"hitX(Y)\t\t";
	    //cout<<"calX(Y)\t\t";
	    //cout<<"hit-track\t";
	    //cout<<"cal-track\t";
	    //cout<<"calY\t\t";
	    //cout<<"hitZ\t";
	    cout<<"chi2\t";
	    cout<<endl;
	    for(Int_t layer=0;layer<NumHitLayer;layer++)
	      {
		cout<<trk1->GetHitID(layer);
		//cout<<"\t"<<layer;
		cout<<"\t"<<setw(8)<<trk1->GetDriftLength(layer);
		if(trk1->GetPosition(0)>-9000)
		  {
		    cout<<"\tX:"<<setw(8)<<trk1->GetPosition(0);
		    //cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer);
		    //cout<<"\t"<<setw(8)<<trk1->GetPosition(0)+trk1->GetAngle(0)*trk1->GetHitZPosition(layer);
		    //cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer)-trk1->GetPosition(0);
		    //cout<<"\t"<<setw(8)<<trk1->GetAngle(0)*trk1->GetHitZPosition(layer);
		  }
		else
		  {
		    cout<<"\tY:"<<setw(8)<<trk1->GetPosition(1);
		    //cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer);
		    //cout<<"\t"<<setw(8)<<trk1->GetPosition(1)+trk1->GetAngle(1)*trk1->GetHitZPosition(layer);
		    //cout<<"\t"<<setw(8)<<trk1->GetHitXPosition(layer)-trk1->GetPosition(1);
		    //cout<<"\t"<<setw(8)<<trk1->GetAngle(1)*trk1->GetHitZPosition(layer);
		  }
	      
		//cout<<"\t"<<setw(8)<<trk1->GetAngle(0);
		//cout<<"\t"<<setw(5)<<trk1->GetHitZPosition(layer);
		cout<<"\t"<<setw(8)<<trk1->GetChi2();

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

    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1);
    c1->cd(1);
    tree->Draw("abs(bdc1calX-bdc1trackX):bdc1tdc>>h1(120,1640,1760,200,0,10)","bdc1layer==0&&bdc1isX==1","colz");
    c1->cd(2);
    tree->Draw("abs(bdc1calY-bdc1trackY):bdc1tdc>>h2(120,1640,1760,200,0,10)","bdc1layer==2&&bdc1isX==0","colz");

    TCanvas *c2 = new TCanvas("c2","c2",1200,600);
    c2->Divide(2,1);
    c2->cd(1);
    tree->Draw("abs(bdc1hitX-bdc1trackX):bdc1tdc>>h3(120,1640,1760,200,0,10)","bdc1layer==0&&bdc1isX==1","colz");
    c2->cd(2);
    tree->Draw("abs(bdc1hitY-bdc1trackY):bdc1tdc>>h4(120,1640,1760,200,0,10)","bdc1layer==2&&bdc1isX==0","colz");

    TCanvas *c3 = new TCanvas("c3","c3",1200,600);
    c3->Divide(2,1);
    c3->cd(1);
    tree->Draw("bdc1calX-bdc1hitX>>h5(100,-10,10)","bdc1layer==0&&bdc1isX==1");
    gPad->SetLogy();
    c3->cd(2);
    tree->Draw("bdc1calY-bdc1hitY>>h6(100,-10,10)","bdc1layer==2&&bdc1isX==0");
    gPad->SetLogy();
  
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
