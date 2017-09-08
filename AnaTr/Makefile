TARGET = AnaTree

#OBJS = AnaTree.o AnaMINOS.o AnaNEBULA.o AnaBeamPla.o AnaBDC.o AnaHOD.o AnaFDC.o AnaNeuLAND.o AnaCoin.o
OBJS = AnaTree.o AnaBeamPla.o AnaBDC.o AnaHOD.o AnaFDC.o AnaNeut.o AnaCoin.o AnaNEBULA.o AnaNeuLAND.o AnaDALI.o AnaCATANA.o AnaPPAC.o Julien.o AnaBrho.o AnaSAMURAI.o MatCalc.o AnaBeamPID.o AnaFragPID.o AnaRelE.o loadPar.o AnaNEBslew.o

CODESOURCE =  $(TARTSYS)

ROOTCFLAGS  = $(shell root-config --cflags)
ROOTLIBS    = $(shell root-config --libs) -lGX11 -lXMLParser
CLIBS = -L$(CODESOURCE)/lib -lanacore  -lanaroot -lananadeko -lanabrips -lanacatana -lanadali -lanasamurai -lanacore
CFLAGS =  -O3  -funroll-loops  -fgcse-lm  -fgcse-sm  -fschedule-insns -Wall
INCLUDES = -I$(CODESOURCE)/include 

CC = g++-4.9

#.SUFFIXES: .C .o

all: $(TARGET) 

$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(ROOTCFLAGS) $(INCLUDES) $(CLIBS) $(ROOTLIBS)

%.o: %.C 
	$(CC) -c $^ -o $@ $(CFLAGS) $(ROOTCFLAGS) $(INCLUDES)

clean:
	rm -f $(TARGET) $(OBJS)

install:
	ln -sf $(PWD)/AnaTree ../../AnaTree 
#AnaTree.o: AnaTree.C AnaTree.H
#AnaMINOS.o: AnaMINOS.C AnaMINOS.H
#AnaNEBULA.o: AnaNEBULA.C AnaNEBULA.H
#AnaBeamPla.o: AnaBeamPla.C AnaBeamPla.H
#AnaBDC.o: AnaBDC.C AnaBDC.H
