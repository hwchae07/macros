#ifndef JULIEN_H
#define JULIEN_H

class Julien {

public:
  Julien(){}
  ~Julien(){}


  void ReadParameters(std::string var, std::ifstream &is, std::vector<Double_t> &res);
  void ReadParameters(std::string var, std::ifstream &is, std::vector<Int_t> &res);
  void ReadParameters(std::string var, std::ifstream &is, Double_t &res);
  void ReadParameters(std::string var, std::ifstream &is, Int_t &res);
  void PrintParameters(std::string var, std::vector<Double_t> res);
  void PrintParameters(std::string var, std::vector<Int_t> res);
  Bool_t ReadParameterMultiDimFit();
  Double_t GetRigidityFitResult(Double_t *x);
  Double_t GetPathLengthFitResult(Double_t *x);

private:
  // rigidity
  Int_t      r_gNVariables    ;
  Int_t      r_gNCoefficients ;
  Double_t   r_gDMean         ;
  std::vector<Double_t>  r_gXMin;
  std::vector<Double_t>  r_gXMax;
  std::vector<Double_t>  r_gCoefficient;
  std::vector<Int_t   >  r_gPower;

  // path length
  Int_t      l_gNVariables    ;
  Int_t      l_gNCoefficients ;
  Double_t   l_gDMean         ;
  std::vector<Double_t>  l_gXMin;
  std::vector<Double_t>  l_gXMax;
  std::vector<Double_t>  l_gCoefficient;
  std::vector<Int_t   >  l_gPower;
};

#endif
