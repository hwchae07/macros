void GetHist(const char* filenames, Int_t idStart = 1, Int_t idMax = 144){

  TChain* tree = new TChain("nebula","nebula");
  tree->Add(filenames);

  TFile* file = new TFile("root/nebula_hpc.root","RECREATE");

  for (Int_t id = idStart ; id <= idMax ; id++){
    for (Int_t comb1 = 0 ; comb1 < 4 ; comb1++){
      for (Int_t comb2 = 0 ; comb2 < 4 ; comb2++){
	Int_t hpc1, hpc2;
	if (id < 121) { // for NEBULA
	  hpc1 = comb1 +     (id - 1)/60 * 8;
	  hpc2 = comb2 + 4 + (id - 1)/60 * 8;}
	else { // for VETO
	  hpc1 = comb1 +     (id - 121)/12 * 8;
	  hpc2 = comb2 + 4 + (id - 121)/12 * 8;}
	tree->Draw(Form("nebulaTD - nebulaTU>>h%d_%d%d(100,-30,30)",id,comb1,comb2),
		   Form("nebulaID == %d && hpcID[%d] && hpcID[%d]",id,hpc1,hpc2),"GOFF");
	cout << "ID: " << id << ", COMB.: (" << hpc1 << ", " << hpc2 << "). Saved." << endl;}}}
  
  file->Write();
  file->Close();
  delete tree;}

void FitHist(const char* filename, Int_t idStart = 1, Int_t idMax = 144, Double_t fitRange = 1, Bool_t flagMon = true){

  TFile* file = new TFile(filename);

  TCanvas* can;
  if (flagMon) {
    can = new TCanvas("FitHist","FitHist",800,600);
    can->Divide(4,4);}

  Double_t fitRes[144][4][4];
  for (Int_t id = 0 ; id < 144 ; id++)
    for (Int_t comb1 = 0 ; comb1 < 4 ; comb1++)
      for (Int_t comb2 = 0 ; comb2 < 4 ; comb2++)
	fitRes[id][comb1][comb2] = 0;

  for (Int_t id = idStart ; id <= idMax ; id++){
    if (flagMon) can->Clear("D");
    for (Int_t comb1 = 0 ; comb1 < 4 ; comb1++){
      for (Int_t comb2 = 0 ; comb2 < 4 ; comb2++){
	Int_t hpc1, hpc2;
	if (id < 121) { // for NEBULA
	  hpc1 = comb1 +     (id - 1)/60 * 8;
	  hpc2 = comb2 + 4 + (id - 1)/60 * 8;}
	else { // for VETO
	  hpc1 = comb1 +     (id - 121)/12 * 8;
	  hpc2 = comb2 + 4 + (id - 121)/12 * 8;}
	TH1 *h1;
	file->GetObject(Form("h%d_%d%d",id,comb1,comb2),h1);
	h1->SetTitle(Form("ID: %d, Comb:(%d,%d);T_D - T_U [ns];",id,comb1,comb2));
	if (flagMon) can->cd(4*comb1 + comb2 + 1);
	Double_t maxX = h1->GetBinCenter(h1->GetMaximumBin());
	if (flagMon) h1->Fit("gaus","Q","",maxX-fitRange,maxX+fitRange);
	else  h1->Fit("gaus","Q0","",maxX-fitRange,maxX+fitRange);
	TF1 *f1;
	f1 = h1->GetFunction("gaus");
	if (f1) {
	  fitRes[id-1][comb1][comb2] = f1->GetParameter(1); // mean value
	  cout << "ID: " << id << ", COMB.: (" << hpc1 << ", " << hpc2 << "). Fitted." << endl;}
	else {
	  fitRes[id-1][comb1][comb2] = 0; // mean value
	  cout << "ID: " << id << ", COMB.: (" << hpc1 << ", " << hpc2 << "). No Data." << endl;}
	if (flagMon) can->Update();}}}
  
  // result to txt
  ofstream fout("nebula_hpc_fit.txt");
  for (Int_t id = idStart ; id <= idMax ; id++){
    for (Int_t comb1 = 0 ; comb1 < 4 ; comb1++){
      for (Int_t comb2 = 0 ; comb2 < 4 ; comb2++)
	fout << fitRes[id-1][comb1][comb2] << " ";}
    fout << endl;}
  fout.close();
  //  file->Close();
}

void FitGraph(const char* filename, Int_t idStart = 1, Int_t idMax = 144, Bool_t flagMon = true){

  Int_t nOmitComb = 3;
  Int_t omitComb[3] = { 10, 22, 32 };
  
  Double_t hpcY[16] = {870, 470, 82.5, -327.5, 
		       600, 200, -205, -915,
		       870, 470, 77.5, -330, 
		       597.5, 197.5, -202.5, -905};
  Double_t hpcZ[4] = { 13555.2, 14365.2, 14401.2, 15211.2 };
  Double_t nebulaY[16];
  Double_t nebulaZ[4] = { 13895.2, 14025.2, 14741.2, 14871.2 };
  Double_t vetoZ[4] = { 13620.2, 13635.2, 14466.2, 14482.2 };


  Double_t fitRes[144][4][4];
  for (Int_t i = 0 ; i < 144 ; i++)
    for (Int_t j = 0 ; j < 4 ; j++)
      for (Int_t k = 0 ; k < 4 ; k++)
	fitRes[i][j][k] = 0;

  ifstream fin("nebula_hpc_fit.txt");
  Int_t id = 0;
  TString line;
  while( line.ReadLine(fin) ){
    TString token;
    Int_t from = 0;
    Int_t j = 0;
    while( line.Tokenize(token,from) ){
      fitRes[id][j/4][j%4] = token.Atof(); j++;}
    id++;}
  fin.close();

  TCanvas *can;
  if (flagMon) {
    can = new TCanvas("FitGraph","FitGraph",800,600);
    can->Divide(6,5);}

  Int_t nPoint;
  Double_t x[16], y[16];
  Double_t result[144][2];
  for (Int_t id = idStart ; id <= idMax ; id++){   
    if (flagMon && id == 121) {
      can->Clear(); can->Divide(6,4);}
    for (Int_t comb1 = 0 ; comb1 < 4 ; comb1++){
      for (Int_t comb2 = 0 ; comb2 < 4 ; comb2++){
	Int_t hpc1, hpc2;
	if (id < 121) { // for NEBULA
	  hpc1 = comb1 +     (id - 1)/60 * 8;
	  hpc2 = comb2 + 4 + (id - 1)/60 * 8;
	  nebulaY[comb1*4+comb2] = hpcY[hpc1] + (hpcY[hpc2]-hpcY[hpc1])/
	    (hpcZ[hpc2 / 4]-hpcZ[hpc1 / 4]) * (nebulaZ[(id - 1)/30]-hpcZ[hpc1 / 4]);}
	else { // for VETO
	  hpc1 = comb1 +     (id - 121)/12 * 8;
	  hpc2 = comb2 + 4 + (id - 121)/12 * 8;
	  nebulaY[comb1*4+comb2] = hpcY[hpc1] + (hpcY[hpc2]-hpcY[hpc1])/
	    (hpcZ[(id - 121)/12 * 2 + 1]-hpcZ[(id - 121)/12 * 2]) 
	    * (vetoZ[(id - 121)/12 * 2 + (id - 121) % 2]-hpcZ[(id - 121)/12 * 2]);}}}

    nPoint = 0;
    for (Int_t i = 0 ; i < 16 ; i++){
      if (fitRes[id-1][i/4][i%4] == 0) continue; // no fitting data
      Bool_t omitNow = false;
      for (Int_t j = 0 ; j < nOmitComb ; j++)
	if ((omitComb[j] / 10 * 4 + omitComb[j] % 10) == i) omitNow = true;
      if (omitNow) continue; // omit this combination
      x[nPoint] = fitRes[id-1][i/4][i%4];
      y[nPoint] = nebulaY[i];
      nPoint++;}
    
    TGraph *gr = new TGraph(nPoint,x,y);
    gr->SetTitle(Form("ID: %d;dT [ns];Y [mm]",id));
    if (flagMon){
      if ( id < 121 ) can->cd((id-1) % 30 + 1);
      else can->cd((id-121) + 1);
      gr->Fit("pol1","Q","");
      gr->Draw("A*");
      can->Update();}
    else
      gr->Fit("pol1","Q0","");
    for (Int_t k = 0 ; k < 2 ; k++) {
      TF1* f1 = gr->GetFunction("pol1");
      result[id-1][k] = f1->GetParameter(k);}}

  // result to txt
  ofstream fout("nebula_hpc_res.txt");
  for (Int_t id = idStart ; id <= idMax ; id++){
    for (Int_t k = 0 ; k < 2 ; k++) {
      fout << result[id-1][k] << " ";}
    fout << endl;}
  fout.close();
}
