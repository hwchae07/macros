{
  //PID cut//
  TFile *pidcut = new TFile("./cut/pid_32ne_physics.root","r");
  TCutG *b34na = (TCutG*)pidcut->Get("b34na");
  TCutG *f32ne = (TCutG*)pidcut->Get("f32ne");
  f32ne->SetVarX("fragAoZ");
  f32ne->SetVarY("fragZ");
  //PID cut//

  TH1 *hist1 = new TH1D("hist1","hist1",100,0,5);
  TH1 *hsum1 = new TH1D("hsum1","hsum1",100,0,5);
  
  for(Int_t runNum = 275;runNum<=294;runNum++)
    //for(Int_t runNum = 295;runNum<=301;runNum++)
    {
      TFile *file = new TFile(Form("./root/run%04d.root.BeamPla",runNum),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("EventNum");
      tree->AddFriend("SAMURAI",Form("./root/run%04d.root.SAMURAI",runNum));
      tree->GetFriend("SAMURAI")->BuildIndex("EventNum");
      tree->AddFriend("FDC",Form("./root/run%04d.root.FDC",runNum));
      tree->GetFriend("FDC")->BuildIndex("EventNum");
      tree->AddFriend("BDC",Form("./root/run%04d.root.BDC",runNum));
      tree->GetFriend("BDC")->BuildIndex("EventNum");
      tree->AddFriend("HOD",Form("./root/run%04d.root.HOD",runNum));
      tree->GetFriend("HOD")->BuildIndex("EventNum");
      tree->AddFriend("Coin",Form("./root/run%04d.root.Coin",runNum));
      tree->GetFriend("Coin")->BuildIndex("EventNum");
      tree->AddFriend("Neut",Form("./root/run%04d.root.Neut",runNum));
      tree->GetFriend("Neut")->BuildIndex("EventNum");
  
      //Beam PID//
      Double_t tof713_offset = 213.73055 +352.345;
      Double_t fl7_13 = (3957.1684+3957.2192) / 2.;    //cm
      Double_t brho0 = 7.8412;    //Tm
      tree->SetAlias("tofc713",Form("TOF7_13+%lf",tof713_offset));
      tree->SetAlias("betaF7",Form("%lf/tofc713/29.9792458",fl7_13));
      tree->SetAlias("betaF5","betaF7*1.005");
      tree->SetAlias("betaF13","betaF7*0.9");
      tree->SetAlias("dEFac_beam","TMath::Log(4866 * betaF13 * betaF13) - TMath::Log(1-betaF13*betaF13) - betaF13*betaF13");
      tree->SetAlias("beamZ","1+10*TMath::Sqrt(icbEloss/dEFac_beam)*betaF13");
      tree->SetAlias("beamMomU"," 931.494*betaF5/TMath::Sqrt(1-betaF5*betaF5)"); // mom/c = gamma*mass*beta [MeV/u/c]
      tree->SetAlias("beamAoZ",Form("(1+f5X/3300)*%lf * 299.792458 / beamMomU ",brho0));
      //Beam PID//
  
      //Target size//
      Double_t bdcWidth = 120.;
      Double_t bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
      Double_t bdc1Tgt = bdc2Tgt + 1000;
      Double_t tgtZ = -4978.89+488;
      tree->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
      tree->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
      TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";  
      //Target size//

      //HOD time offset//
      /*ifstream file1("./dat/hod_offset_brho_scan.dat");
      Int_t dummy[25];
      Double_t offset[25];
      for(Int_t i=0;i<=24;i++)
      file1>>dummy[i]>>offset[i];*/
      TString tof_hod = "";
      tof_hod += "hodTA-T13";
      /*tof_hod += "+";
      tof_hod += offset[0];
      for(Int_t id=1;id<=24;id++)
	{
	  tof_hod += "+";
	  tof_hod += offset[id];
	  tof_hod += Form("*(hodID==%d)",id);
	  }*/
      tree->SetAlias("tof_hod",tof_hod);
      //HOD time offset//

      //FDC z position//
      Double_t fdc1Z = -2888.82;
      Double_t alpha = -59.92/180.*TMath::Pi();
      Double_t fdc2Z = 4164.51;    //3726.51+438;
      tree->SetAlias("fdc2zz",Form("%lf*cos(%lf) - (fdc2X+715.01)*sin(%lf)",fdc2Z,alpha,alpha));
      tree->SetAlias("fdc2xx",Form("%lf*sin(%lf) + (fdc2X+715.01)*cos(%lf)",fdc2Z,alpha,alpha));
      //FDC z position//

      //intersection//
      tree->SetAlias("a1",Form("tan( fdc1A+%lf )",0.));
      tree->SetAlias("b1",Form("fdc1X - %lf*a1",fdc1Z));
      tree->SetAlias("a2",Form("tan(fdc2A+%lf)",alpha));
      tree->SetAlias("b2","fdc2xx - fdc2zz*a2");
      tree->SetAlias("interZ","(b2-b1)/(a1-a2)");
      tree->SetAlias("interX","(a1*b2-b1*a2)/(a1-a2)");
      //intersection//

      //hod position//
      Double_t hodZ = 4999.17;
      tree->SetAlias("hodX",Form("fdc2X + 715.01 + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
      tree->SetAlias("hodzz",Form("%lf*cos(%lf) - (hodX)*sin(%lf)",hodZ,alpha,alpha));
      tree->SetAlias("hodxx",Form("%lf*sin(%lf) + (hodX)*cos(%lf)",hodZ,alpha,alpha));
      tree->SetAlias("theta",Form("abs(%lf + fdc1A - fdc2A)",alpha));
      //hod position//

      //fragment PID//
      tree->SetAlias("FL1",Form("sqrt( (interX-tgtX)**2 + (interZ-%lf)**2 )",tgtZ));
      tree->SetAlias("FL2","sqrt( (interX-hodxx)**2 + (interZ-hodzz)**2 )");
      tree->SetAlias("FL","FL1 + FL2 + 1000*brho/2.9*theta - 2*1000*brho/2.9*tan(theta/2)");

      tree->SetAlias("beta","FL/tof_hod/299.792458");
      tree->SetAlias("dEFac_frag","TMath::Log(15795.98 * beta * beta) - TMath::Log(1-beta*beta) - beta*beta");
      tree->SetAlias("fragZ","TMath::Sqrt(hodQA/TMath::Sqrt((1+TMath::Tan(fdc2A)*TMath::Tan(fdc2A))*(1+fdc2B*fdc2B))/dEFac_frag)*beta");
      tree->SetAlias("fragMomUfromHOD"," 931.494*beta/TMath::Sqrt(1-beta*beta)"); // mom/c = gamma*mass*beta [MeV/u/c]
      tree->SetAlias("fragAoZ","brho * 299.792458 / fragMomUfromHOD ");
      //fragment PID//

      //coincidence with neutron//
      tree->SetAlias("trig_neut","coinTrigger==3||coinTrigger==5||coinTrigger==7");
      //coincidence with neutron//
      
  
      //Neutron//
      Double_t m_neut = massMeV(1,0);
      Double_t amu_neut = massAMU(1,0);
      tree->SetAlias("FL_neut","sqrt((neutX-tgtX)**2+(neutY-tgtY)**2+neutZ**2)");
      tree->SetAlias("neutBeta","FL_neut/(neutTA-T13)/299.792458");
      tree->SetAlias("neutGamma","1/sqrt(1-neutBeta**2)");
      tree->SetAlias("neutE",Form("neutGamma*%lf",m_neut));
      tree->SetAlias("neutKE",Form("(neutGamma-1)*%lf/%lf",m_neut,amu_neut));
      tree->SetAlias("neutP",Form("neutGamma*%lf*neutBeta",m_neut));
      tree->SetAlias("neutPx","neutP*(neutX-tgtX)/FL_neut");
      tree->SetAlias("neutPy","neutP*(neutY-tgtY)/FL_neut");
      tree->SetAlias("neutPz","neutP*(neutZ)/FL_neut");
      //Neutron//
  
      //32Ne 4momentum//
      Double_t m32ne = massMeV(32,10);
      Double_t amu32ne = massAMU(32,10);
      Double_t zz = 10;
      tree->SetAlias("fragP",Form("brho*299.792458*%lf",zz));
      tree->SetAlias("fragPx","fragP*fdc1ex");
      tree->SetAlias("fragPy","fragP*fdc1ey");
      tree->SetAlias("fragPz","fragP*fdc1ez");
      tree->SetAlias("fragE",Form("sqrt(fragP**2 + %lf**2)",m32ne));
      tree->SetAlias("fragKE",Form("(fragE-%lf)/%lf",m32ne,amu32ne));
      //32Ne 4momentum//


      tree->SetAlias("relE",Form("sqrt((fragE+neutE)**2 - (fragPx+neutPx)**2 -(fragPy+neutPy)**2 - (fragPz+neutPz)**2) - %lf - %lf",m_neut,m32ne));


      
      tree->Draw("relE>>hist1(100,0,5)","b34na&&f32ne&&trig_neut&&neutVETO==0"&&cut_reaction);
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
      
      
    }
  Double_t Ntot = 379201.;
  Double_t areaDensity = 1.08*TMath::Power(10,23) * TMath::Power(10,-27);    //mb

  Double_t binWidth = hsum1->GetBinWidth(1);
  hsum1->Scale(1/Ntot/areaDensity/binWidth);
  
  hsum1->Draw();
  
}
