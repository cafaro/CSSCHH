CC=clang++
# optimization: on mac we use -Os on linux and other unix systems -O3
CFLAGS=-Os -std=c++11
LDFLAGS=-lm
SOURCES=csschh.cc spacesaving.cc prng.cc rand48.cc massdal.cc
EXECUTABLE=csschh

all:
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(SOURCES) $(LDFLAGS)

clean:
	rm -f *~ $(EXECUTABLE)
