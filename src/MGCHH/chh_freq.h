#include "frequent.h"

//using namespace std;

#define min(x,y)	fmin(x, y)
#define max(x,y)	fmax(x, y)


typedef struct input_parameters {
	unsigned long stream_size;						// N
	unsigned int primary_distinct_items;   	//Universe U1
	unsigned int secondary_distinct_items; 	//Universe U2
	unsigned int s1, s2; 	// Sketch dimensions	
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
	unsigned long long int item;
	unsigned int count;
} ChhCounter;

typedef struct {
	std::map <unsigned long long int,ChhCounter> primary;
	std::map <unsigned long long int,ChhCounter> correlated;
} frequent_items;

void 	frequent_Print				(frequent_items *fi);
    
    
