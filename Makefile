CXX = g++

SOURCES = det_track.cxx kinema.h kinema.cxx mkdata.h mkdata.cxx \
	anatpc.h anatpc.cxx para.h kinema_doi.cpp kinema_doi.hpp database.cpp\
	database.hpp nuclear.cpp nuclear.hpp dataset.cpp

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --glibs)
DEBUG = -Wno-unused-but-set-variable -Wno-write-strings

CFLAGS = -c -O3 -lm $(ROOTFLAGS)
CFLAGS += -std=c++0x

CFLAGS2  = -O3 -lm
CFLAGS2 += -std=c++0x

all: det_track #add

.cpp.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<
.cxx.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<
.c.o:
	$(CXX) $(CFLAGS) $(DEBUG) -I$(GARFIELD_HOME)/Include $<

det_track: det_track.o database.o kinema.o mkdata.o anatpc.o kinema_doi.o database.o nuclear.o dataset.o
	$(CXX) $(CFLAGS2) $(DEBUG) det_track.o mkdata.o kinema.o anatpc.o kinema_doi.o database.o nuclear.o dataset.o \
	-o det_track \
	$(GARFIELD_HOME)/Library/libGarfield.a \
	-lm -lgfortran $(ROOTLIBS)
#	rm -f *.o

add: add.o
	$(CXX) $(CFLAGS) $(ROOTLIBS) $(DEBUG) -o add add.o

database.o: database.hpp
kinema.o: kinema.h para.h
mkdata.o: mkdata.h para.h
anatpc.o: anatpc.h para.h
kinema_doi.o: kinema_doi.hpp
database.o: database.hpp
nuclear.o: nuclear.hpp

clean:
	rm -f det_track add *.o
