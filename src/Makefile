﻿CXX     = g++
CXXFLAGS = -Wall -g -O2 -std=c++0x -I../include
LDLIBS = -lpthread

JELLYFISH = /home/mitra/tools/jellyfish-2.1.4
JELLYFISH_CFLAGS = -I$(JELLYFISH)/include -I$(JELLYFISH)
JELLYFISH_LIBS = -L$(JELLYFISH)/.libs -Wl,-Bstatic -ljellyfish-2.0 -Wl,-Bdynamic
CXXFLAGS += $(JELLYFISH_CFLAGS)
LDLIBS += $(JELLYFISH_LIBS)

OBJS =  perfectread_main.o

all: perfectread.bin

perfectread: $(OBJS)
	$(CXX) $(OBJS) -o perfectread.bin $(LDLIBS)

perfectread_main.o: perfectread_main.cc perfectread.hpp  perfectread_main_cmdline.hpp  qual_mer_dna.hpp  sequence_handler.hpp
	$(CXX) $(CXXFLAGS) -c perfectread_main.cc

clean:
	rm -f *.o perfectread.bin

