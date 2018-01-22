{
  gROOT->ProcessLine("Double_t *center = new Double_t[2]");
  gROOT->ProcessLine("Double_t center1[24],center2[24]");
  gROOT->ProcessLine(".L ./macros/brho_scan_hodQ.cc");

  for(Int_t id = 1;id<=23;id++)
    {
      gROOT->ProcessLine(Form("center = brho_scan_hodQ(%d)",id));
      gROOT->ProcessLine(Form("center1[%d] = center[0]",id));
      gROOT->ProcessLine(Form("center2[%d] = center[1]",id));
    }

  gROOT->ProcessLine("for(Int_t id=1 ; id<24 ; id++){");
  gROOT->ProcessLine("cout<<center1[id]<<endl<<center2[id]<<endl;}");
}
