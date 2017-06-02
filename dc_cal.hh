
TArtSAMURAIParameters *fBDCParameters;
TArtCalibBDC1Hit *fCalibBDC1Hit;
TArtCalibBDC2Hit *fCalibBDC2Hit;
TArtCalibBDC1Track *fCalibBDC1Track;
TArtCalibBDC2Track *fCalibBDC2Track;

Int_t bdc1idNum;
Int_t bdc1id;
/*
Int_t bdc1layer[100];
Int_t bdc1tdc[100];
Double_t bdc1dl[100];
*/

std::vector<int> bdc1layer;
std::vector<int> bdc1tdc;
std::vector<double> bdc1dl;

/*
Int_t *bdc1layer;
Int_t *bdc1tdc;
Double_t *bdc1dl;
*/
void LoadTDCDistribution(const char *filename);
void analysis(Int_t number=1);
void getDL(Int_t number=1);
void LoadTDCDistribution(const char *filename);
