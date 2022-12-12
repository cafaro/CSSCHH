CC=clang++
# optimization: on mac we use -Os on linux and other unix systems -O3
CFLAGS=-Os -std=c++11
LDFLAGS=-lm
SOURCES=chh_freq.cpp frequent.cpp prng.cpp rand48.cpp massdal.cc
EXECUTABLE=chh_freq

all:
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(SOURCES) $(LDFLAGS)

clean:
	rm -f *~ $(EXECUTABLE)
