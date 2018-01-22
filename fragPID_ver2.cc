{

  
  TFile *fcut = new TFile("./cut/beamCut.root","r");
  TCutG *b32na = (TCutG*)fcut->Get("b32na");
  TCutG *b33na = (TCutG*)fcut->Get("b33na");
  TCutG *b34na = (TCutG*)fcut->Get("b34na");
  TCutG *b30ne = (TCutG*)fcut->Get("b30ne");
  TCutG *b31ne = (TCutG*)fcut->Get("b31ne");
  TCutG *b32ne = (TCutG*)fcut->Get("b32ne");
  TCutG *b29f = (TCutG*)fcut->Get("b29f");
  fcut->Close();
  
  // Int_t runNum = 275;
  // TString filename = Form("run%04d.root",runNum);

  // TFile *file = new TFile(Form("./root/%s.BeamPID",filename.Data()),"r");
  // TTree *tree = (TTree*)file->Get("BeamPID");
  // tree->BuildIndex("EventNum");
  // tree->AddFriend("BDC",Form("./root/%s.BDC",filename.Data()));
  // tree->GetFriend("BDC")->BuildIndex("EventNum");
  // tree->AddFriend("FDC",Form("./root/%s.FDC",filename.Data()));
  // tree->GetFriend("FDC")->BuildIndex("EventNum");
  // tree->AddFriend("HOD",Form("./root/%s.HOD",filename.Data()));
  // tree->GetFriend("HOD")->BuildIndex("EventNum");
  // tree->AddFriend("SAMURAI",Form("./root/%s.SAMURAI",filename.Data()));
  // tree->GetFriend("SAMURAI")->BuildIndex("EventNum");


  // //Target size//                                                                                                             
  // Double_t bdcWidth = 120.;
  // Double_t bdc2Tgt = 5366 - 4479 - bdcWidth/2.;
  // Double_t bdc1Tgt = bdc2Tgt + 1000;
  // Double_t tgtZ = -4978.89+488;

  // tree->SetAlias("tgtX",Form("(bdc2X - bdc1X)/1000*%lf + bdc1X",bdc1Tgt));
  // tree->SetAlias("tgtY",Form("(bdc2Y - bdc1Y)/1000*%lf + bdc1Y",bdc1Tgt));
  // TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";
  // //Target size//                                                                                                             

  // TString tof_hod = "";
  // tof_hod += "hodTA-T13";
  // tree->SetAlias("tof_hod",tof_hod);

  // //FDC z position//                                                                                                          
  // Double_t fdc1Z = -2888.82;
  // Double_t alpha = -59.92/180.*TMath::Pi();
  // Double_t fdc2Z = 4164.51;    //3726.51+438;                                                                                 

  // tree->SetAlias("fdc2zz",Form("%lf*cos(%lf) - (fdc2X+715.01)*sin(%lf)",fdc2Z,alpha,alpha));
  // tree->SetAlias("fdc2xx",Form("%lf*sin(%lf) + (fdc2X+715.01)*cos(%lf)",fdc2Z,alpha,alpha));
  // //FDC z position//                                                                                                          

  // //intersection//                                                                                                            
  // tree->SetAlias("a1",Form("tan( fdc1A+%lf )",0.));
  // tree->SetAlias("b1",Form("fdc1X - %lf*a1",fdc1Z));
  // tree->SetAlias("a2",Form("tan(fdc2A+%lf)",alpha));
  // tree->SetAlias("b2","fdc2xx - fdc2zz*a2");
  // tree->SetAlias("interZ","(b2-b1)/(a1-a2)");
  // tree->SetAlias("interX","(a1*b2-b1*a2)/(a1-a2)");
  // //intersection//                                                                                                            


  // Double_t hodZ = 4999.17;
  // tree->SetAlias("hodX",Form("fdc2X + 715.01 + (%lf-%lf)*tan(fdc2A)",hodZ,fdc2Z));
  // //tree->SetAlias("hodX","2667.14-2667.14/24/2-hodID*2667.14/24");                                                           
  // tree->SetAlias("hodzz",Form("%lf*cos(%lf) - (hodX)*sin(%lf)",hodZ,alpha,alpha));
  // tree->SetAlias("hodxx",Form("%lf*sin(%lf) + (hodX)*cos(%lf)",hodZ,alpha,alpha));
  // tree->SetAlias("theta",Form("abs(%lf + fdc1A - fdc2A)",alpha));

  // tree->SetAlias("FL1",Form("sqrt( (interX-tgtX)**2 + (interZ-%lf)**2 )",tgtZ));
  // tree->SetAlias("FL2","sqrt( (interX-hodxx)**2 + (interZ-hodzz)**2 )");
  // tree->SetAlias("FL","FL1 + FL2 + 1000*brho/2.9*theta - 2*1000*brho/2.9*tan(theta/2)");

  // tree->SetAlias("beta","FL/tof_hod/299.792458");



  // tree->SetAlias("dEFac","TMath::Log(15795.98 * beta * beta) - TMath::Log(1-beta*beta) - beta*beta");

  // tree->SetAlias("fragZ","TMath::Sqrt(hodQ/TMath::Sqrt( (1+tan(fdc2A)**2 + tan(fdc2B)**2) )/dEFac)*beta");
  // /*  
  //     if(align==true)
  //     tree->SetAlias("fragZ","TMath::Sqrt(hodQ/TMath::Sqrt( (1+tan(fdc2A)**2 + tan(fdc2B)**2) )/dEFac)*beta");
  //     else
  //     tree->SetAlias("fragZ","TMath::Sqrt(hodQA/TMath::Sqrt( (1+tan(fdc2A)**2 + tan(fdc2B)**2) )/dEFac)*beta");
  // */	

  // tree->SetAlias("fragMomUfromHOD"," 931.494*beta/TMath::Sqrt(1-beta*beta)"); // mom/c = gamma*mass*beta [MeV/u/c]            
  // tree->SetAlias("fragAoZ","brho * 299.792458 / fragMomUfromHOD ");


  // //tree->Draw("fragZ:fragAoZ>>histPID(200,2.0,3.5,200,3.5,5)",cut_34na&&cut_reaction&&"hodNum==1&&trig_neut","colz");

  // TCanvas *c1 = new TCanvas("c1","c1",1200,600);
  // c1->Divide(2,1);
  // c1->cd(1);
  // tree->Draw("fragZ:fragAoZ>>hist1(400,2.0,3.5,400,3.5,5.0)","b34na","colz");
      
  
  
  //c1->cd(2);
  TChain *chain = new TChain("FragPID","FragPID");
  //chain->Add("./root/run0275.root.FragPID");
  chain->Add("./root/run027[5-9].root.FragPID");
  chain->Add("./root/run028[0-9].root.FragPID");
  chain->Add("./root/run029[0-4].root.FragPID");
  
  TCut cut_reaction = "TMath::Sqrt( TMath::Power(tgtX,2) + TMath::Power(tgtY,2)) < 40";
  TCut trig_neut = "coinTrigger==3||coinTrigger==5||coinTrigger==7";
  TCut cut_f32ne = "abs(fragAoZ-3.34)<0.06&&abs(fragZ-3.95)<0.15";

  Double_t f32ne_x[9] = {3.32112,3.27802,3.26455,3.3292,3.35614,3.39116,3.40194,3.35075,3.32112};
  Double_t f32ne_y[9] = {4.04911,4.01403,3.90242,3.84503,3.83227,3.89286,4.00446,4.05548,4.04911};

  TCutG *f32ne = new TCutG("f32ne",9,f32ne_x,f32ne_y);
  f32ne->SetVarX("fragAoZ");
  f32ne->SetVarY("fragZ");

  
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2);
  c1->cd(1);
  chain->Draw("fragZ:fragAoZ>>h1(400,2.0,3.5,400,3.5,5.0)","b34na"&&cut_reaction,"colz");
  h1->SetTitle("Fragment PID for{}^{34}Na beam w/o neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
  f32ne->Draw("SAME");
  c1->cd(2);
  chain->Draw("fragZ:fragAoZ>>h2(200,2.0,3.5,200,3.5,5.0)","b34na"&&cut_reaction&&trig_neut,"colz");
  h2->SetTitle("Fragment PID for{}^{34}Na beam w/ neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
  h2->GetZaxis()->SetRangeUser(1,100);

  
  c1->cd(3);
  chain->Draw("fragZ:fragAoZ>>h3(400,2.0,3.5,400,3.5,5.0)","b34na"&&cut_reaction&&"f32ne","colz");
  h3->SetTitle("Fragment PID for{}^{34}Na beam w/o neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
  c1->cd(4);
  chain->Draw("fragZ:fragAoZ>>h4(200,2.0,3.5,200,3.5,5.0)","b34na"&&cut_reaction&&trig_neut&&"f32ne","colz");
  h4->SetTitle("Fragment PID for{}^{34}Na beam w/ neutron coincidence;A/Z (arbitrary);Z (arbitrary)");
  h4->GetZaxis()->SetRangeUser(1,100);
  
  
  /*
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    c1->Divide(2,2);
    c1->cd(1);
    chain->Draw("fragZ:fragAoZ>>h1(400,2.0,3.5,400,3.5,5.0)","b34na","colz");
    c1->cd(2);
    chain->Draw("fragZ:fragAoZ>>h2(400,2.0,3.5,400,3.5,5.0)","b33na","colz");
    c1->cd(3);
    chain->Draw("fragZ:fragAoZ>>h3(400,2.0,3.5,400,3.5,5.0)","b32ne","colz");
    c1->cd(4);
    chain->Draw("fragZ:fragAoZ>>h4(400,2.0,3.5,400,3.5,5.0)","b31ne","colz");
  */
}
