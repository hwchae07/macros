#ifndef MatCalc_H
#define MatCalc_H
#include "TMatrixD.h"

class MatCalc{

public:
  MatCalc(const char *filename, Double_t x0, Double_t a0);
  ~MatCalc(){}

  void LoadMatrix(const char* filename);
  Double_t CalcDelta(TMatrixD& input,TMatrixD& output, Int_t xora);
  Double_t Brho;
  
private:
  TMatrixD zero;
  TMatrixD first;
  TMatrixD second[4];

  Double_t x_offset;
  Double_t a_offset;

};

#endif
