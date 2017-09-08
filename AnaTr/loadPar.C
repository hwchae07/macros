#include <iostream>
#include <fstream>
#include "TString.h"
#include "loadPar.H"

/*
void test()
{
  double **tdc1 = loadParU(6);
  double **tdc2 = loadParD(6);
  cout<<tdc1[0][0]<<" "<<tdc1[0][1]<<endl;
  cout<<tdc1[1][0]<<" "<<tdc1[1][1]<<endl;

  cout<<tdc2[0][0]<<" "<<tdc2[0][1]<<endl;
  cout<<tdc2[1][0]<<" "<<tdc2[1][1]<<endl;
  
}
*/

double **loadParU(int order)
{
  TString numberingList[6] = {"1st","2nd","3rd","4th","5th","6th"};
  TString numbering = numberingList[order-1];

  std::ifstream fin(Form("../../dat/NEBULA/par_%s_time_vs_channel_u.dat",numbering.Data()));

  int idNum = 144;
  //double tdcPar[144][2];

  
  double **tdcPar = new double*[idNum];
  for(int id=0;id<idNum;id++)
    tdcPar[id] = new double[order+1];
  
  for(int id=0;id<144;id++)
    {
      for(int i=0;i<=order;i++)
	fin>>tdcPar[id][i];
    }

  fin.close();
  
  return tdcPar;
  
  for(int id=0;id<idNum;id++)
    delete [] tdcPar[id];

  delete [] tdcPar;
  	
}

double **loadParD(int order)
{
  TString numberingList[6] = {"1st","2nd","3rd","4th","5th","6th"};
  TString numbering = numberingList[order-1];

  std::ifstream fin(Form("../../dat/NEBULA/par_%s_time_vs_channel_d.dat",numbering.Data()));

  int idNum = 144;
  //double tdcPar[144][2];

  
  double **tdcPar = new double*[idNum];
  for(int id=0;id<idNum;id++)
    tdcPar[id] = new double[order+1];
  
  for(int id=0;id<144;id++)
    {
      for(int i=0;i<=order;i++)
	fin>>tdcPar[id][i];
    }


  fin.close();

  return tdcPar;

  for(int id=0;id<idNum;id++)
    delete [] tdcPar[id];

  delete [] tdcPar;
  	
}
