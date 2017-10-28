CC = g++
CFLAGS = -c

CXXFLAGS = -DLINUXVERS -g -Wall -I$(ROOTSYS)/include -DWITH_DEBUG -I$(ROOTSYS) -I$(ROOTSYS)/test/
LDFLAGS =
ROOTLIBS = $(shell root-config --glibs) 
LIBS = -lnsl -lcrypt -ldl $(ROOTLIBS) -lMinuit $(XLIB) -lXt -lX11 -lXpm  

SOURCES = root_anna_pl.cpp
HEADER = root_anna_pl.h


# #*********************************************
# # First method of makefile

OBJECTS = $(SOURCES:.c=.o) $(ROOTSYS)/test/libEvent.so
EXECUTABLE = main

all: $(SOURCES) $(EXECUTABLE)

#	$(CC) -lm -o $(EXECUTABLE) $(OBJECTS) $@
$(EXECUTABLE): $(OBJECTS)
	$(CC) -lm -o $(EXECUTABLE) $(OBJECTS) $(LIBS) $(CXXFLAGS)

#	$(CC) -lm $(EXECUTABLE) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o  $(EXECUTABLE) $(OBJECTS) $(LIBS) $(HEADER) $(CXXFLAGS) $@ 

clean:
	rm -f $(EXECUTABLE) *.o


