#ifndef ANATEST_H
#define ANATEST_H

#include <csignal>

#include "TFile.h"

// Variable for signal handling
bool flagEnd;

// Signal Handling Procedure
void SigHandler(int param){
  flagEnd = true;
}

void PrintUsage();
void PrintTitle();
int SelectModule();

TFile* LoadRootFile(const char* name, const char* detName, bool readOnly = true);

#endif
