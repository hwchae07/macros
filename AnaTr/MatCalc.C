#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TPRegexp.h"
#include "TMatrixD.h"

#include "MatCalc.H"

void test(const char* filename)
{
  MatCalc *load = new MatCalc(filename,0,0);
}

MatCalc::MatCalc(const char* filename, Double_t x0, Double_t a0)
{
  x_offset = x0;
  a_offset = a0;
  LoadMatrix(filename);
}

void MatCalc::LoadMatrix(const char* filename){
  std::ifstream fin(filename);
  if(!fin.is_open())
    return;

  zero.ResizeTo(6,1);    //determine matrix size? 6 by 1 matrix? vector
  first.ResizeTo(6,6);   //6 by 6 matrix?
  for(Int_t i = 0 ; i<4 ; i++)    //second order matrix... why 4?
    second[i].ResizeTo(6,6);

  TString line;
  Int_t lineNum=0;

  while(!(line.ReadLine(fin)).eof()){
    if(line[0]=='c'||line[0]=='C')
      continue;
    lineNum++;

    if(lineNum==1)
      Brho = line.Atof();
    else if(lineNum==2)
      zero(0,0) = line.Atof();
    else if(lineNum==3)
      zero(1,0) = line.Atof();
    else if(lineNum <=9)
      {
	TStringToken token(line," ");    //string is sequentially placed in token.
	for(Int_t i=0 ; token.NextToken() ; i++)
	  first(lineNum-4,i) = token.Atof();
      }
    else
      {
	Int_t mat_index = (lineNum - 10) / 6;
	Int_t fir_index = (lineNum - 10) % 6;
	TStringToken token(line," ");
	for(Int_t i=0 ; token.NextToken() ; i++)
	  second[mat_index](fir_index,i) = token.Atof();
      }

    //std::cout<<lineNum<<" "<<line.Atof()<<endl;
  }

  zero(0,0) += x_offset;
  zero(1,0) += a_offset;

  /*
  std::cout<<zero(0,0)<<std::endl;
  std::cout<<zero(1,0)<<std::endl;
  for(Int_t i=0;i<6;i++)
    {
      std::cout<<std::endl;
      for(Int_t j=0; j<6;j++)
	std::cout<<first(i,j)<<"   ";
    }
  */
  fin.close();
  
}

Double_t MatCalc::CalcDelta(TMatrixD& input, TMatrixD& output, Int_t xora)
{
  
  Double_t coeff[3];
  Double_t delta=0;
  
  TMatrixD calc_out(6,1);
  calc_out.Mult(first,input);
  calc_out += zero;

  //coeff[0] = first(0,0) * input(0,0) + first(0,1) * input(1,0) + zero(0,0) - output(0,0);
  coeff[0] = calc_out(xora,0) - output(xora,0);
  coeff[1] = first(xora,5);
  coeff[2] = 0;


  delta = -coeff[0] / coeff[1];

 
  
  //2nd order//
  for(Int_t i=0 ; i<4 ; i++)
    {
      for(Int_t j=0 ; j<6 ; j++)
	{
	  for(Int_t k=j ; k<6 ; k++)
	    {
	      calc_out(i,0) += (second[i](j,k) * input(j,0) * input(k,0));
	    }
	}
    }

  coeff[0] = calc_out(xora,0) - output(xora,0);
  coeff[1] = first(xora,5);
  coeff[2] = second[xora](5,5);
  
  for(Int_t i=0 ; i<5 ; i++)
  coeff[1] += (second[xora](i,5) * input(i,0));
  coeff[2] = second[xora](5,5);      
  //2nd order//
  
  Double_t D = coeff[1]*coeff[1] - 4 * coeff[0] * coeff[2];
  if(D > 0)
    delta = (-coeff[1] + TMath::Sqrt(D))/(2*coeff[2]);
  else
    delta = (-coeff[0] / coeff[1]);

  return delta;
}


