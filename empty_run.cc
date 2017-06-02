{

  //gROOT->Reset();

  //Int_t runNum = 296;

  TH2 *hist1 = new TH2D("hist1","hist1",400,2.8,3.5,400,3,5);
  TH2 *hsum1 = new TH2D("hsum1","hsum1",400,2.8,3.5,400,3,5);
  
  for(Int_t runNum=296;runNum<=301;runNum++)
    {
      TString filename = Form("run%04d.root",runNum);
  
      TFile *file = new TFile(Form("./root/%s.BeamPla",filename.Data()),"r");
      TTree *tree = (TTree*)file->Get("BeamPla");
      tree->BuildIndex("EventNum");
      tree->AddFriend("BDC",Form("./root/%s.BDC",filename.Data()));
      tree->GetFriend("BDC")->BuildIndex("EventNum");
      tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
      tree->GetFriend("FDC")->BuildIndex("EventNum");
      tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
      tree->GetFriend("HOD")->BuildIndex("EventNum");
      //tree->AddFriend("Brho",Form("./root/%s.Brho",filename.Data()));
      //tree->GetFriend("Brho")->BuildIndex("EventNum");
      tree->AddFriend("SAMURAI",Form("./root/%s.SAMURAI",filename.Data()));
      tree->GetFriend("SAMURAI")->BuildIndex("EventNum");
      tree->AddFriend("Coin",Form("./root/%s.Coin",filename.Data()));
      tree->GetFriend("Coin")->BuildIndex("EventNum");


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
	file1>>dummy[i]>>offset[i];
      */
      TString tof_hod = "";
      tof_hod += "hodTA-T13";
      /*tof_hod += "+";
      tof_hod += offset[0];

      for(Int_t id=1;id<=24;id++)
	{
	  tof_hod += "+";
	  tof_hod += offset[id];
	  tof_hod += Form("*(hodID==%d)",id);
	}
      */
      tree->SetAlias("tof_hod",tof_hod);
      //HOD time offset//

      //FDC z position//
      Double_t fdc1Z = -2888.82;
      Double_t alpha = -59.92/180.*TMath::Pi();
      Double_t fdc2Z = 3726.51;
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

      //HOD position//
      Double_t hodZ = 5004.17;
      tree->SetAlias("hodX",Form("fdc2X + 715.01 + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));    //HOD x position by extarpolation
      tree->SetAlias("hodzz",Form("%lf*cos(%lf) - (hodX)*sin(%lf)",hodZ,alpha,alpha));    //absolute x position of HOD
      tree->SetAlias("hodxx",Form("%lf*sin(%lf) + (hodX)*cos(%lf)",hodZ,alpha,alpha));    //absolute z position of HOD
      tree->SetAlias("theta",Form("abs(%lf + fdc1A - fdc2A)",alpha));
      //HOD position//

      //Flight length of target to HOD//
      tree->SetAlias("FL1",Form("sqrt( (interX-tgtX)**2 + (interZ-%lf)**2 )",tgtZ));
      tree->SetAlias("FL2","sqrt( (interX-hodxx)**2 + (interZ-hodzz)**2 )");
      tree->SetAlias("FL","FL1 + FL2 + 1000*brho/2.9*theta - 2*1000*brho/2.9*tan(theta/2)");
      //Flight length of target to HOD//

      //beta reconstruction//
      tree->SetAlias("betaHOD","FL/tof_hod/299.792458");
      //tree->SetAlias("betaHOD","9495.06/tof_hod/299.792458");  // 9495.06 = 4490.89 + 5004.17
      tree->SetAlias("gammaHOD","1/sqrt(1-betaHOD**2)");
      //beta reconstruction//

      //fragment PID//
      tree->SetAlias("dEFacFrag","TMath::Log(15795.98 * betaHOD * betaHOD) - TMath::Log(1-betaHOD*betaHOD) - betaHOD*betaHOD");
      tree->SetAlias("fragZ","TMath::Sqrt(hodQA/TMath::Sqrt((1+TMath::Tan(fdc2A)*TMath::Tan(fdc2A))*(1+fdc2B*fdc2B))/dEFacFrag)*betaHOD");
  
      tree->SetAlias("fragMomUfromHOD"," 931.494*betaHOD/TMath::Sqrt(1-betaHOD*betaHOD)"); // mom/c = gamma*mass*beta [MeV/u/c]
      tree->SetAlias("fragAoZ","brho * 299.792458 / fragMomUfromHOD ");
      //fragment PID//

      //fragment energy from Brho for only 32Ne//
      tree->SetAlias("mom_brho","brho*299.792458*10");
      tree->SetAlias("energy_brho","(sqrt(mom_brho**2 +29839.6979**2 )-29839.6979)/32");
      //fragment energy from Brho for only 32Ne//


  
  
  
      tree->SetAlias("tofc","TOF3_13-3*dT5");
      tree->SetAlias("elossc","icbEloss-0.1*tofc-40");
   
      TCut cut_na = "abs(elossc-22.5)<2.5";
      TCut cut_ne = "abs(elossc-17.8)<2.2";
      TCut cut_32ne = cut_ne&&"abs(tofc+379)<3";
      TCut cut_34na = cut_na&&"abs(tofc+388)<3";
  
      TCut cut_central = "abs(f5X)<10";


      //TOF F7~F13 offset setting//
      //TCanvas *c1 = new TCanvas("c1","c1",600,600);
      tree->Draw("TOF7_13>>h1(100,-354,-351)",cut_central&&cut_32ne,"goff");
      TH1 *h1;
      gDirectory->GetObject("h1",h1);
      TF1 *fit1 = new TF1("fit1","gaus",-354,-351);
      h1->Fit(fit1,"rq","",-353,-352);
      //c1->Close();
  
      Double_t offset_pla = 213.73055 - fit1->GetParameter(1);
      tree->SetAlias("tofc713",Form("TOF7_13+%lf",offset_pla));
      Double_t par[4] = {6.84494091796875000000E+03,
			 -7.73498382568359375000E+01,
			 3.06769594550132751465E-01,
			 -4.17593779275193810463E-04};
      //TOF F7~F13 offset setting//
  
      //TOF F7~F13 to Energy//
      TString energy = "";
      energy += Form("%lf",par[0]);
      for(Int_t i=1;i<4;i++)
	energy += Form("+ (%lf) * TMath::Power(tofc713,%d)",par[i],i);

      Double_t fl7_13 = (3957.1684+3957.2192) / 2. ;    //cm
      Double_t brho0 = 7.8412;    //Tm

      tree->SetAlias("energy",energy);
      //tree->SetAlias("gamma","energy/931.494+1");
      //tree->SetAlias("beta","sqrt(1-1/gamma**2)");
      //TOF F7~F13 to Energy//
  
      tree->SetAlias("beta",Form("%lf/tofc713/29.9792458",fl7_13));
      tree->SetAlias("gamma","1./sqrt(1-beta**2)");
      tree->SetAlias("betaF5","beta*1.005");
      tree->SetAlias("betaF13","beta*0.9");
  
      tree->SetAlias("dEFac","TMath::Log(4866 * betaF13 * betaF13) - TMath::Log(1-betaF13*betaF13) - betaF13*betaF13");
      tree->SetAlias("beamZ","1+10*TMath::Sqrt(icbEloss/dEFac)*betaF13");

  
      tree->SetAlias("beamMom","29839.6979*beta/TMath::Sqrt(1-beta*beta)"); // mom/c = gamma*mass*beta [MeV/u/c]
      tree->SetAlias("beamMomU"," 931.494*betaF5/TMath::Sqrt(1-betaF5*betaF5)"); // mom/c = gamma*mass*beta [MeV/u/c]
      tree->SetAlias("beamAoZ",Form("(1+f5X/3300)*%lf * 299.792458 / beamMomU ",brho0));

      TCut cut_beam = "abs(beamZ-10)<2 && abs(beamAoZ-3)<0.5";

      TFile *cutfile = new TFile("./cut/pid_empty_run.root","r");
      TCutG *b33na,*b34na,*b30ne,*b31ne,*b32ne,*b29f;
      b34na = (TCutG*)cutfile->Get("b34na");
      b33na = (TCutG*)cutfile->Get("b33na");
      b30ne = (TCutG*)cutfile->Get("b30ne");
      b31ne = (TCutG*)cutfile->Get("b31ne");
      b32ne = (TCutG*)cutfile->Get("b32ne");
      b29f = (TCutG*)cutfile->Get("b29f");

  
      /*  
	  TCanvas *c2 = new TCanvas("c2","c2",600,800);
	  c2->Divide(1,2);
	  c2->cd(1);
	  tree->Draw("betaHOD>>beta1(100,0.6,0.7)",cut_32ne);
	  c2->cd(2);
	  tree->Draw("beta>>beta2(100,0.6,0.7)",cut_32ne);
      */

      /*
	TCanvas *c3 = new TCanvas("c3","c3",1500,1000);
	c3->Divide(3,2);
	c3->cd(1);
	tree->Draw("fragZ:fragAoZ>>fpid1(200,3,4,200,3.5,5.5)","b33na","colz");
	c3->Update();
	c3->cd(2);
	tree->Draw("fragZ:fragAoZ>>fpid2(200,3,4,200,3.5,5.5)","b34na","colz");
	c3->Update();
	c3->cd(3);
	tree->Draw("fragZ:fragAoZ>>fpid3(200,3,4,200,3.5,5.5)","b29f","colz");
	c3->Update();
	c3->cd(4);
	tree->Draw("fragZ:fragAoZ>>fpid4(200,3,4,200,3.5,5.5)","b30ne","colz");
	c3->Update();
	c3->cd(5);
	tree->Draw("fragZ:fragAoZ>>fpid5(200,3,4,200,3.5,5.5)","b31ne","colz");
	c3->Update();
	c3->cd(6);
	tree->Draw("fragZ:fragAoZ>>fpid6(200,3,4,200,3.5,5.5)","b32ne","colz");
	c3->Update();
      */
      /*
	TCanvas *c4 = new TCanvas("c4","c4",1200,600);
	c4->Divide(2,1);
	c4->cd(1);
	tree->Draw("brho*299.792458/(gamma*beta*931.494):fdc2X>>pid11(200,200,500,200,2,4)","b32ne","colz");
	c4->cd(2);
	tree->Draw("brho*299.792458/(gammaHOD*betaHOD*931.494):fdc2X>>pid12(200,200,500,200,2,4)","b32ne","colz");
      */
      /*
	TCanvas *c5 = new TCanvas("c5","c5",1000,500);
	c5->Divide(2,1);
	c5->cd(1);
	tree->Draw("fragZ:brho*299.792458/(gamma*beta*931.494)>>pid1(200,2.8,3.6,200,3.5,5.5)","","colz");
	c5->cd(2);
	tree->Draw("fragZ:brho*299.792458/(gammaHOD*betaHOD*931.494)>>pid2(200,2.8,3.6,200,3.5,5.5)","","colz");
      */

      //coincidence with neutron//
      tree->SetAlias("trig_neut","coinTrigger==3||coinTrigger==5||coinTrigger==7");
      //coincidence with neutron//
      
      
      tree->Draw("fragZ:fragAoZ>>hist1(400,2.8,3.5,400,3,5)","b34na&&trig_neut","goff");
      gDirectory->GetObject("hist1",hist1);
      hsum1->Add(hist1);
  
    }

  hsum1->Draw("colz");
}
