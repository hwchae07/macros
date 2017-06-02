#include <iostream>

#include "AnaMINOS.H"

#include "TArtTrackMINOSData.hh"

#include "TTree.h"
#include "TXMLNode.h"
#include "TDOMParser.h"

AnaMINOS::AnaMINOS():
  AnaModule("MINOS"),
  fMINOSParameters(NULL),
  fCalibMINOS(NULL),
  fAnalyzedMINOS(NULL),
  fTrackMINOS(NULL),
  fVertexMINOS(NULL),
  minosXv( -100 ),
  minosYv( -100 ),
  minosZv( -100 ),
  minosTrX( -100 ),
  minosTrA( -100 ),
  minosTrY( -100 ),
  minosTrB( -100 ),
  MINOSthresh( 25.  ),
  TimeBinElec( 20.  ),
  VDrift    (   4.0701),
  Tshaping  ( 333.  ),
  DelayTrig (1200.0 ),
  Tgt_Length( 152.0 ),
  Beta      (   0.6 ),
  Delta_Beta(   0.0 ),
  Pos_offset( -7.0 ),
  DALIOffset(  0.0 ) {
  ;}

AnaMINOS::~AnaMINOS(){    
  DeleteAll();}

void AnaMINOS::InitParameter(){
  fMINOSParameters = TArtMINOSParameters::Instance();
  fMINOSParameters->LoadParameters((char*)"db/MINOS.xml");
  LoadParameters("db/MINOS_Analysis.xml");
  parLoaded = true;}

void AnaMINOS::InitDetector(){
  if (!parLoaded) return;
  fCalibMINOS = new TArtCalibMINOS;
  fAnalyzedMINOS = new TArtAnalyzedMINOS(fCalibMINOS);
  fTrackMINOS = new TArtTrackMINOS;
  fVertexMINOS = new TArtVertexMINOS;

  fAnalyzedMINOS->SetConfig(VDrift, TimeBinElec, DelayTrig);
  detLoaded = true;}

void AnaMINOS::Analysis(){
  anaFlag = true;
  if (!detLoaded) return;

  fCalibMINOS->ClearData();
  fAnalyzedMINOS->ClearData();
  fTrackMINOS->ClearData();
  fVertexMINOS->ClearData();

  fCalibMINOS->ReconstructData();
  fAnalyzedMINOS->ReconstructData();
 if (fAnalyzedMINOS->GetNumAnalyzedMINOS() <= 10) {
    anaFlag = false; return; }
  fTrackMINOS->ReconstructData();
  if (fTrackMINOS->GetTrackNumMINOS() < 1) {
    anaFlag = false; return; }
  fVertexMINOS->ReconstructData();

  minosTrack = fTrackMINOS->GetTrackNumMINOS();
  minosXv = fVertexMINOS->GetXv();
  minosYv = fVertexMINOS->GetYv();
  minosZv = fVertexMINOS->GetZv();

  TArtTrackMINOSData* tr = fTrackMINOS->GetTrackMINOS(0);  
  minosTrX = tr->GetPar_x0();
  minosTrA = tr->GetPar_Ax();
  minosTrY = tr->GetPar_y0();
  minosTrB = tr->GetPar_Ay();
}

void AnaMINOS::DeleteAll(){
  //  if (fMINOSParameters) delete fMINOSParameters;
  if (fCalibMINOS)    delete fCalibMINOS;
  if (fAnalyzedMINOS) delete fAnalyzedMINOS;
  if (fTrackMINOS)    delete fTrackMINOS;
  if (fVertexMINOS)   delete fVertexMINOS;
}

void AnaMINOS::SetTree(){ 
  AnaModule::SetTree(); 
  tree->Branch("minosTrack",&minosTrack,"minosTrack/I");
  tree->Branch("minosXv",&minosXv,"minosXv/D");
  tree->Branch("minosYv",&minosYv,"minosYv/D");
  tree->Branch("minosZv",&minosZv,"minosZv/D");
  tree->Branch("minosTrX",&minosTrX,"minosTrX/D");
  tree->Branch("minosTrA",&minosTrA,"minosTrA/D");
  tree->Branch("minosTrY",&minosTrY,"minosTrY/D");
  tree->Branch("minosTrB",&minosTrB,"minosTrB/D");}

bool AnaMINOS::LoadParameters(const char *xmlfile)
{
  TDOMParser domParser;
  domParser.SetValidate(false);
  Int_t parsecode = domParser.ParseFile(xmlfile);
  if(parsecode < 0){
    std::cerr << domParser.GetParseCodeMessage(parsecode) << std::endl;
    return false;
  }
  TXMLNode* node = domParser.GetXMLDocument()->GetRootNode();

  if(strcmp(node->GetNodeName(), "dataroot")) return false;
  node = node->GetChildren();
  for(; node; node = node->GetNextNode()){
    if(node->GetNodeType() != TXMLNode::kXMLElementNode) continue; // Element Node    
    ParseParaList(node->GetChildren());
  }
  return true;}

void AnaMINOS::ParseParaList(TXMLNode *node){
  for(; node; node = node->GetNextNode()){
    if(node->GetNodeType() != TXMLNode::kXMLElementNode) continue; // Element Node
    if(strcmp(node->GetNodeName(), "MINOSthresh") == 0){
      MINOSthresh = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "TimeBinElec") == 0){
      TimeBinElec = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "VDrift") == 0){
      VDrift = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "Tshaping") == 0){
      Tshaping = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "DelayTrig") == 0){
      DelayTrig = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "Tgt_Length") == 0){
      Tgt_Length = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "Beta") == 0){
      Beta = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "Delta_Beta") == 0){
      Delta_Beta = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "Pos_offset") == 0){
      Delta_Beta = atof(node->GetText());
    } else if(strcmp(node->GetNodeName(), "DALIOffset") == 0){
      DALIOffset = atof(node->GetText());
    }
  }
}
