void NebulaTree()
{
  TArtSAMURAIParameters *setup = new TArtSAMURAIParameters();
  setup->LoadParameter("./db/NEBULA.xml");
  setup->LoadParameter("./db/NEBULAHPC.xml");

  TArtEventStore *estore = new TArtEventStore();

  TString ridffile = "./sdaq14/nebula0027.ridf";

  estore->Open(ridffile.Data());

  TArtCalibNEBULA *calibnebula = new TArtCalibNEBULA;
  TArtCalibNEBULAHPC *calibhpc = new TArtCalibNEBULAHPC;

  Int_t neve = 0;

  TTree *tree = new TTree("nebula","nebula");

  /*
  nebulaID = new Int_t[144];
  nebulaX = new Double_t[144];
  nebulaY = new Double_t[144];
  nebulaZ = new Double_t[144];
  */

  Int_t hpcNum = 0;
  Int_t hpcID[145]={0,};
  Int_t nebulaNum = 0;
  Int_t nebulaID[145]={0,};
  Double_t nebulaX[145]={0,};
  Double_t nebulaY[145]={0,};
  Double_t nebulaZ[145]={0,};

  tree->Branch("nebulaID",nebulaID,"nebulaID[144]/I");
  tree->Branch("nebulaX",nebulaX,"nebulaX[144]/D");
  tree->Branch("nebulaY",nebulaY,"nebulaY[144]/D");
  tree->Branch("nebulaZ",nebulaZ,"nebulaZ[144]/D");

  tree->Branch("hpcID",hpcID,"hpcID[144]/I");

  while(estore->GetNextEvent()&&neve<10000)
    {
      if (neve%1000==0) std::cout<<"\r"<<"event : "<<neve<<std::flush;

      calibnebula->ClearData();
      calibnebula->ReconstructData();

      calibhpc->ClearData();
      calibhpc->ReconstructData();

      hpcNum = 0;
      nebulaNum = 0;

      for(Int_t i=0 ; i<calibnebula->GetNumNEBULAPla() ; i++)
	{
	  TArtNEBULAPla *pla = calibnebula->GetNEBULAPla(i);
	  if(pla->GetHit())
	    {
	      //cout<<"ID,X : "<<pla->GetID()<<","<<pla->GetPos(0)<<endl;
	      nebulaID[nebulaNum] = pla->GetID();
	      nebulaX[nebulaNum] = pla->GetPos(0);
	      nebulaY[nebulaNum] = pla->GetPos(1);
	      nebulaZ[nebulaNum] = pla->GetPos(2);

	      nebulaNum++;
	    }
	}

      for(Int_t i=0 ; i<calibhpc->GetNumNEBULAHPC() ; i++)
	{
	  TArtNEBULAHPC *hpc = calibhpc->GetNEBULAHPC(i);
	  hpcID[hpcNum] = hpc->GetID();

	  hpcNum++;
	}
      tree->Fill();
      ++neve;

      /*
      Int_t NumNEBULAPla = calibnebula->GetNumNEBULAPla();
      for(Int_t i=0 ; i<NumNEBULAPla ; i++)
	{
	  TArtNEBULAPla *pla = calibnebula->GetNEBULAPla(i);

	  if(pla->GetHit())
	    {
	      nebulaID[i] = pla->GetID();
	      nebulaX[i] = pla->GetPos(0);
	      nebulaY[i] = pla->GetPos(1);
	      nebulaZ[i] = pla->GetPos(2);
	    }
	}
      tree->Fill();
      ++neve;
      */

    }

  estore->ClearData();


  /*
  delete nebulaID;
  delete nebulaX;
  delete nebulaY;
  delete nebulaZ;
  */
}
