#include "spacesaving.h"

//using namespace std;

#define min(x,y)	fmin(x, y)
#define max(x,y)	fmax(x, y)


typedef struct input_parameters {
	unsigned long stream_size;						// N
	unsigned int primary_distinct_items;   	//Universe U1
	unsigned int secondary_distinct_items; 	//Universe U2
	unsigned int k1, k2; 	// Sketch dimensions
	double zipfpar, hwzpar;
	double phi1, eps1, phi2, eps2;
	unsigned int seed;
	int hseed;
	bool verbose;
	bool csvFormat;
	char* testFile;
	bool withTestFile;
} input_parameters;

typedef struct {
	std::map <SSLItem,SSLCounterdata> primary;
	std::map <SSLItem,SSLCounterdata> correlated;
} frequent_items;

void 	frequent_Print				(frequent_items *fi);


