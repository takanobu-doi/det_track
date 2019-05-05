CXX = g++

SOURCES = det_track.cpp kinema.h kinema.cpp mkdata.h mkdata.cpp \
	anatpc.h anatpc.cpp para.h 12C_3alpha.hpp 12C_3alpha.cpp

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)
DEBUG = -Wno-unused-but-set-variable -Wno-write-strings

CFLAGS = -c -O3 -lm $(ROOTFLAGS)
CFLAGS += -std=c++0x

CFLAGS2  = -O3 -lm
CFLAGS2 += -std=c++0x

all: det_track #add

.cc.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<
.cpp.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<
.cxx.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<
.c.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<

det_track: det_track.o kinema.o mkdata.o anatpc.o 12C_3alpha.o
	$(CXX) $(CFLAGS2) $(DEBUG) det_track.o mkdata.o kinema.o anatpc.o 12C_3alpha.o -o det_track \
	$(GARFIELD_HOME)/Library/libGarfield.a \
	-lm -lgfortran $(ROOTLIBS)
#	rm -f *.o

add: add.o
	$(CXX) $(CFLAGS) $(ROOTLIBS) $(DEBUG) -o add add.o

det_track.o: para.h
kinema.o: kinema.h para.h
mkdata.o: mkdata.h para.h
anatpc.o: anatpc.h para.h
12C_3alpha.o: 12C_3alpha.hpp para.h

clean:
	rm -f det_track add *.o
clean~:
	$(RM) *~
