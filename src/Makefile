CXX     = g++
CXXFLAGS = -Wall -g -O2 -std=c++0x -I../include
LDLIBS = -lpthread

JELLYFISH = /home/mitra/tools/jellyfish-2.1.4
JELLYFISH_CFLAGS = -I$(JELLYFISH)/include -I$(JELLYFISH)
JELLYFISH_LIBS = -L$(JELLYFISH)/.libs -Wl,-Bstatic -ljellyfish-2.0 -Wl,-Bdynamic
CXXFLAGS += $(JELLYFISH_CFLAGS)
LDLIBS += $(JELLYFISH_LIBS)

OBJS =  perfectread_main.o

all: perfectread

perfectread: $(OBJS)
	$(CXX) $(OBJS) -o perfectread $(LDLIBS)

clean:
	rm -f *.o perfectread

