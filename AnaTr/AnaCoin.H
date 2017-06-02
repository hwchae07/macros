#ifndef ANACOIN_H
#define ANACOIN_H

#include "AnaModule.H"

#include "TArtCalibCoin.hh"

class AnaCoin : public AnaModule {

public:

  AnaCoin();
  virtual ~AnaCoin();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:
    TArtCalibCoin* fCalibCoin;

  // for trees
  Int_t coinTrigger;
};

#endif
