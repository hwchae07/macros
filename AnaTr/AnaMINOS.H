#ifndef ANAMINOS_H
#define ANAMINOS_H

#include "AnaModule.H"

#include "TXMLNode.h"

#include "TArtMINOSParameters.hh"
#include "TArtCalibMINOS.hh"
#include "TArtAnalyzedMINOS.hh"
#include "TArtTrackMINOS.hh"
#include "TArtVertexMINOS.hh"

class AnaMINOS : public AnaModule {

public:

  AnaMINOS();
  virtual ~AnaMINOS();

  void InitParameter();
  void InitDetector();
  void Analysis();
  
  void SetTree();

  void DeleteAll();

private:

  bool LoadParameters(const char*);
  void ParseParaList(TXMLNode*);
  
  TArtMINOSParameters* fMINOSParameters;
  
  TArtCalibMINOS* fCalibMINOS;
  TArtAnalyzedMINOS* fAnalyzedMINOS;
  TArtTrackMINOS* fTrackMINOS;
  TArtVertexMINOS* fVertexMINOS;

  // for trees
  Int_t minosTrack;
  Double_t minosXv, minosYv, minosZv;
  Double_t minosTrX, minosTrA, minosTrY, minosTrB;

  // for paramters
  Double_t  MINOSthresh;
  Double_t  TimeBinElec;
  Double_t  VDrift;
  Double_t  Tshaping;
  Double_t  DelayTrig;
  Double_t  Tgt_Length; // target length
  Double_t  Beta; //at beginning of target
  Double_t  Delta_Beta;// beta loss in target 
  Double_t  Pos_offset;//Offset between target & TPC
  Double_t  DALIOffset;
};

#endif
