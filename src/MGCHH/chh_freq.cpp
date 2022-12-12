/********************************************************************
 Approximate correlated frequent items in a data stream
 
 *********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "massdal.h"
#include "prng.h"
#include "chh_freq.h"

using namespace std;

#define LOOPCOUNT 1

typedef struct  {
	int chh_num_primary_candidates; //number of frequent primary candidates
	int chh_num_primary_exact;      //number of true frequent primary items
	int chh_num_correlated_candidates; //number of frequent correlated candidates
	int chh_num_correlated_exact;      //number of true frequent correlated items
	int exact_num_primary; 			//number of exact frequent primary items
	int exact_num_correlated; 		//total number of exact frequent correlated items
	float recall_primary;			//Number of real frequent primary items found / exact number of frequent primary items (recall < 1 implies algorithm fails)
	float precision_primary;		//Number of real frequent primary items found / total frequent primary items found 
	float recall_correlated;		//Number of real frequent correlated items found / exact number of correlated primary items (recall < 1 implies algorithm fails)
	float precision_correlated;		//Number of real frequent correlated items found / total frequent correlated items found 


	float max_rel_error_correlated_candidate;	//Maximum Relative error on frequency estimate for correlated items
	float min_rel_error_correlated_candidate;	//Minimum Relative error (greater than 0) on frequency estimate for correlated items
	float mean_rel_error_correlated_candidate;	//Mean Relative error on frequency estimate for correlated items

	float max_abs_error_correlated_candidate;	//Maximum Relative error on frequency estimate for correlated items
	float mean_abs_error_correlated_candidate;	//Mean Relative error on frequency estimate for correlated items

	float tolerance_eps1; //Experimentally measured eps1 value
	float tolerance_eps2; //Experimentally measured eps2 value

	double update_time;				//Elapsed time for sketch updating
	double query_time;				//Elapsed time for sketch querying
	unsigned int ss_size;		//Number of byte for the sketch data structures
} diagnostic;

void Set_default(input_parameters *p);
void CheckArguments(int argc, char **argv, input_parameters *p) ;
void PrintArguments(input_parameters *p);
unsigned long long int* CreateStream(input_parameters *p);
frequent_items* RunExact(unsigned long long int* s, input_parameters *p);
unsigned long ** AllocateExactCounters(input_parameters *p);
void CheckOutput (input_parameters *p, unsigned long long int *s, frequent_items *e, frequent_items *chh, freq_type* freq, diagnostic *diag);
void PrintDiagnostics(input_parameters *p, diagnostic *diag);
void CHH_Decrement_Randomly (freq_type* freq);
void CHH_Decrement_Least (freq_type* freq);
void CHH_Update(freq_type* freq, unsigned long long int* stream, input_parameters* p);

void frequent_Print(frequent_items *fi);
frequent_items* FindCHH(freq_type *freq, input_parameters *p); 

int main(int argc, char **argv)
{
    diagnostic diag;
    unsigned long long int *stream = NULL;
	input_parameters param;
	frequent_items* exact, *chh;
	freq_type *freq;
		
//	srandom(15);
	
	Set_default(&param);    
    
    CheckArguments(argc,argv, &param);
    
    stream = CreateStream(&param);
	PrintArguments(&param);
        
	printf("Exact computation...\n");
    
    exact = RunExact(stream, &param);

//	printf("Exact Correlated Heavy Hitters\n");
//	frequent_Print(exact);        
    

    diag.update_time = 0;

    for (int j = 0; j < LOOPCOUNT; j++){
	    freq = Freq_Init(param.s1, (void (*)(void*)) &Freq_Destroy );
        
        StartTheClock();


	    CHH_Update(freq, stream, &param);
    	    
		diag.update_time += StopTheClock();
        
        if (j != LOOPCOUNT - 1){
		    Freq_Destroy(freq);
		}
    }

    printf("\nUpdate Time (starttheclock): %f\n\n", diag.update_time);
    
    diag.update_time = diag.update_time/(double)LOOPCOUNT;
	diag.update_time = diag.update_time/(double)param.stream_size;

    StartTheClock();
	chh = FindCHH(freq, &param);
    diag.query_time=StopTheClock();

//	printf("Approximated Correlated Heavy Hitters\n");
//	frequent_Print(chh);        

    CheckOutput(&param, stream, exact, chh, freq, &diag);
	 
    PrintDiagnostics(&param, &diag);
    
    
    delete chh;
    delete exact;
    Freq_Destroy(freq);

    free(stream);

    return 0;
}

void CHH_Update(freq_type* freq, unsigned long long int* stream, input_parameters* p) 
{
	ITEMLIST *item_node, *il;
	GROUP * g;
	freq_type* freq_s;

	for (long i = 0; i < p->stream_size; i++){ 
    	unsigned int x = stream[i] >> 32;
		unsigned int y = stream[i] & 0x00000000FFFFFFFF;
			
		item_node = Freq_Update(freq, x);

		if (item_node) {
				//The primary node has been found (either created or incremented), we can add the secondary item
			if (item_node->value == NULL)
				item_node->value = Freq_Init(p->s2, NULL);
        	Freq_Update((freq_type*)item_node->value, y);
        }
        else {
    	    	//The primary item has not been found neither inserted 
    	    	//we have to remove randomly one of the secondary items from each primary counters

			for (g = freq->groups->nextg; g; g = g->nextg) {
		        il = g->items;
       			if (il) {
		        	do {
		        		freq_s = (freq_type*) il->value;
		        		CHH_Decrement_Randomly(freq_s);
		        		//CHH_Decrement_Least(freq_s);
		            	il = il->nexting;
           			} while (il != g->items);
           		}
		    }
        }
    }
}

void CHH_Decrement_Randomly (freq_type* freq)
{
	int rm;
	
	if (!freq->groups->nextg)
		return;
	rm = random() % (freq->k+1);
	while(freq->items[rm].parentg->diff == 0) {
		rm = (rm+1) % (freq->k+1);
	}
	SubtractCounter(freq->items + rm);
}

void CHH_Decrement_Least (freq_type* freq)
{
	if (freq && freq->groups->nextg && freq->groups->nextg->items)
		Freq_Update(freq, -((int)freq->groups->nextg->items->item));
}

frequent_items* FindCHH(freq_type *freq, input_parameters *p) 
{
	std::vector<FCounter> counters_p, counters_c;
	frequent_items *chh;

	double thresh1, thresh2;
	int i, j;
	
	chh = new frequent_items;
	
    counters_p = Freq_Output(freq);

    thresh1 = (p->phi1 - 1.0/ (double)p->s1) * p->stream_size;
    for (i = 0; i < counters_p.size(); i++) {
   	 	unsigned int x;
   	 	unsigned int f_x;
		x = counters_p[i].item;
   	 	f_x = counters_p[i].count;
   	 	
        if (f_x > thresh1) {
        	ChhCounter xpc;
        	xpc.item = x;
    	    xpc.count = f_x;
			chh->primary[x] = xpc;
			
			counters_c = Freq_Output((freq_type*)counters_p[i].value);
			
			for (j = 0; j < counters_c.size(); j++) {
        		unsigned int y = counters_c[j].item;
		        unsigned int f_xy = counters_c[j].count;
			
				thresh2 = (p->phi2 - 1.0 / (double)p->s2)*(double)f_x - (double)p->stream_size / (double)p->s1;
				if (f_xy > thresh2) {
        			ChhCounter xypc;
        			xypc.item = (unsigned long long int)x << 32 | y;
		    	    xypc.count = f_xy;
				
					chh->correlated[xypc.item] = xypc;
				}
			}
		}
	}

	return chh;
}

void Set_default(input_parameters *p)
{
	double alpha;
	p->stream_size = 1000000;						// N
	p->primary_distinct_items = 1048576;   	//Universe U1
	p->secondary_distinct_items = 1048576;   	//Universe U2
	p->phi1 = 0.01;
	p->phi2 = 0.01;
	p->eps1 = 0.001;
	p->eps2 = 0.001;
	
	alpha = (1.0 + p->phi2)/(p->phi1 - p->eps1);
	if (p->eps1 >= p->eps2 / (2.0 * alpha)){
		p->s1 = ceil ( 2.0 * alpha / p->eps2 );
		p->s2 = ceil( 2.0 / p->eps2 );	
	}
	else {
		p->s1 = ceil (  1.0 / p->eps1 );
		p->s2 = ceil (  1.0 / (p->eps2 - alpha * p->eps1) );
	}
	p->zipfpar = 1.1;
	p->hwzpar = 0;
	p->seed = 16033099;
	p->hseed = 44545;
	p->verbose = false;
	p->csvFormat = false;
	p->testFile = NULL;
	p->withTestFile = false;
	
	return;
}

void CheckArguments(int argc, char **argv, input_parameters *p) 
{    
	double alpha;
    bool failed = false;
    int recomp_bound = 0;
    int recomp_s1s2 = 0;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-ni") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing number of items.\n");
                failed = true;
            }
            p->stream_size = strtol(argv[i], NULL, 10);
        }  else if (strcmp(argv[i], "-tf") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing test file\n");
                failed = true;
            }
            p->testFile = strdup(argv[i]);
            p->withTestFile = true;
        } else if (strcmp(argv[i], "-seed") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing seed.\n");
                failed = true;
            }
            p->seed = (unsigned int)strtol(argv[i], NULL, 10);
        } else if (strcmp(argv[i], "-hseed") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing hash seed.\n");
                failed = true;
            }
            p->hseed = (int) strtol(argv[i], NULL, 10);
            if (p->hseed <= 0)
                p->hseed = static_cast<int>(time(NULL));
        } else if (strcmp(argv[i], "-phi1") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing threshold parameter.\n");
                failed = true;
            }
            p->phi1 = strtod(argv[i], NULL);
            if (p->phi1 <=0 || p->phi1 >= 1) {
            	fprintf(stderr,"phi1 parameters must be a real number in the range (0..1)\n");
            	failed = true;
            } 
        } else if (strcmp(argv[i], "-phi2") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing threshold parameter.\n");
                failed = true;
            }
            p->phi2 = strtod(argv[i], NULL);
            if (p->phi1 <=0 || p->phi1 >= 1) {
            	fprintf(stderr,"phi2 parameters must be a real number in the range (0..1)\n");
            	failed = true;
            } 
        } else if (strcmp(argv[i], "-eps1") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing allowed estimation error.\n");
                failed = true;
            }
            p->eps1 = strtod(argv[i], NULL);
            if (p->eps1 ==0 ) {
            	fprintf(stderr,"eps1 parameters must be a real numer in the range (0..phi1)\n");
            	failed = true;
            }
   			recomp_s1s2 = 1;

        } else if (strcmp(argv[i], "-eps2") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing allowed estimation error.\n");
                failed = true;
            }
            p->eps2 = strtod(argv[i], NULL);
            if (p->eps2 ==0 ) {
            	fprintf(stderr,"eps2 parameters must be a real numer in the range (0..phi2)\n");
            	failed = true;
            }
   			recomp_s1s2 = 1;
        } else if (strcmp(argv[i], "-sk") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing zipfian skew parameter.\n");
                failed = true;
            }
            p->zipfpar = strtod(argv[i], NULL);
        }  else if (strcmp(argv[i], "-hz") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing hurwitz parameter.\n");
                failed = true;
            }
            p->hwzpar = strtod(argv[i], NULL);
        } else if (strcmp(argv[i], "-pdi") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing number of possible secondary distinct items.\n");
                failed = true;
            }
            p->primary_distinct_items = (unsigned int)strtol(argv[i], NULL, 10);
        } else if (strcmp(argv[i], "-sdi") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Missing number of possible secondary distinct items.\n");
                failed = true;
            }
            p->secondary_distinct_items = (unsigned int)strtol(argv[i], NULL, 10);
        } else if (strcmp(argv[i], "-s1") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Number of counters for primary items.\n");
                failed = true;
            }
            p->s1 = (unsigned int)strtol(argv[i], NULL, 10);
            recomp_bound = 1;
            
        } else if (strcmp(argv[i], "-s2") == 0) {
            i++;
            if (i >= argc) {
                fprintf(stderr,"Number of counters for correlated items.\n");
                failed = true;
            }
            p->s2 = (unsigned int)strtol(argv[i], NULL, 10);
            recomp_bound = 1;

        } else if (strcmp(argv[i], "-csv") == 0) {
            p->csvFormat = true;
        } else
            failed = true;
    }

    if (p->zipfpar < 0.0) {PrintArguments(p); failed = true;}
    if (p->phi1 <= 0.0 || p->phi1 >= 1) {PrintArguments(p); failed = true;}
    if (p->phi2 <= 0.0 || p->phi2 >= 1) {PrintArguments(p); failed = true;}
    if (!recomp_bound && (p->eps1 <= 0.0 || p->eps1 >= p->phi1)) {PrintArguments(p); failed = true;}
    if (!recomp_bound && (p->eps2 <= 0.0 || p->eps2 >= p->phi2)) {PrintArguments(p); failed = true;}
    
    
	if (recomp_s1s2) {
		alpha = (1.0 + p->phi2)/(p->phi1 - p->eps1);
		if (p->eps1 >= p->eps2 / (2.0 * alpha)){
			p->s1 = ceil ( 2.0 * alpha / p->eps2 );
			p->s2 = ceil( 2.0 / p->eps2 );	
		}
		else {
			p->s1 = ceil (  1.0 / p->eps1 );
			p->s2 = ceil (  1.0 / (p->eps2 - alpha * p->eps1) );
		}
		if (recomp_bound) 
			printf("\nWarning! options -s1 and -s2 are ignored when used in conjunction with -esp1, -eps2\n\n"); 
	}
	else if (recomp_bound){

		p->eps1 = 1.0 / (float)p->s1;
		p->eps2 = 1.0 / (float)p->s2 + (1.0 + p->phi2)/(p->s1 * p->phi1 - 1.0);			
			
		if (p->eps2 <= 0){
			printf("\x1b[31m Warning ! The size of the first summary s1 does not guarantee a recall of 100%% neither a bound for false positives\x1b[0m \n");
		}
		
		if (p->eps1 >= p->phi1 || p->eps2 >= p->phi2){
			printf("\x1b[31m Warning ! The s1 and s2 values do not guarantee a bound on the false positives \x1b[0m \n"); 	
		}
	}

	if (p->s2 < 1.0/p->phi2){
		printf("\x1b[31m Warning ! The size of the nested summary s2 does not guarantee a recall of 100%% \x1b[0m \n"); 	
	}
    
    if (failed==1)
    {
        printf("\nUsage: %s [operand value] \n\nPossible operands: \n\n",argv[0]);
        printf("-csv = format output as comma separated values. Default = false \n");
        printf("-ni = number of items to process. Default = 100000 \n");
        printf("-sk = parameter of zipf dbn. 0.0 = uniform. 3+ = skewed. Default = 1.1\n");
        printf("-hz = parameter of hurwitz distribution. Default = 0.0\n");
        printf("-phi1 = threshold density for primary heavy hitters.  Default = 0.01\n");
        printf("-phi2 = threshold density for correlated heavy hitters.  Default = 0.01\n");
        printf("-eps1 = max error allowed for primary items estimate (it must be less then the threshold density phi1). Default = 0.001\n");
        printf("-eps2 = max error allowed for correlated items estimate (it must be less then the threshold density phi2). Default = 0.001\n");
        printf("-pdi = number of possible primary distinct items.  Default = 1048575\n");
        printf("-sdi = number of possible secondary distinct items.  Default = 1048575\n");
        printf("-s1 = number of sketch buckets along the primary dimension.\n");
        printf("-s2 = number of sketch buckets along the secondary dimension.\n");
        printf("-seed = seed for input stream generation.  Default = 16033099\n");
        printf("-hseed = seed for hash function generation. Default = 44545\n");
        printf("-tf = file of input dataset. Default = dataset generated\n\n");
        exit(1);
    }
}

void PrintArguments(input_parameters *p)
{
	printf("Input parameters:\n");
	printf("Stream Size: %lu\n",p->stream_size);						// N
	printf("Distinct Items: Primary = %u \t Secondary = %u\n",p->primary_distinct_items, p->secondary_distinct_items);   	//Universe U1
	printf("Thresholds, phi1 = %f \t phi2 = %f\n",p->phi1, p->phi2);
	printf("Error bounds, eps1 = %f \t eps2 = %f\n",p->eps1, p->eps2);
	printf("Sketch size: m1 = %d \t m2 = %d\n",p->s1, p->s2);
	if (!p->withTestFile)
		printf("Input distribution: zipf = %f \t hwz = %f\n", p->zipfpar, p->hwzpar);
	else 
		printf("Real distribution: %s\n",p->testFile);
	printf("Seed = %d \t Hseed = %d\n", p->seed, p->hseed);
	printf("Verbosity: %d\n",p->verbose);
	printf("csvFormat: %d\n",p->csvFormat);
	
}

unsigned long long int* CreateStream(input_parameters *p)
{
	unsigned long long int * stream;
	unsigned long long int x, y;
	
    if (!p->withTestFile) {
        stream = (unsigned long long int *) calloc(p->stream_size, sizeof(unsigned long long int));

        Tools::Random r1 = Tools::Random(p->seed);
        Tools::Random r2 = Tools::Random(p->seed + 1);
        Tools::PRGHurwitz hurwitz1 = Tools::PRGHurwitz(1, p->primary_distinct_items+1, p->zipfpar, p->hwzpar, &r1);
        Tools::PRGHurwitz hurwitz2 = Tools::PRGHurwitz(1, p->secondary_distinct_items+1, p->zipfpar, p->hwzpar, &r2);

        for (long i = 0; i < p->stream_size; i++) {
            x = hurwitz1.nextHurwitzLong();
            y = hurwitz2.nextHurwitzLong();
            stream[i] = (x<<32) | y;
        }
    } else {
		std::ifstream inputFile;

	    inputFile.open(p->testFile, std::ios::binary | std::ios::in);
    	if (!inputFile.is_open()) {
			std::cerr << "Error: bad input file" << std::endl;
			exit (1);
    	}

		inputFile.seekg(0, std::ios::beg);
    	std::streampos begin = inputFile.tellg();
    	inputFile.seekg (0, std::ios::end);
    	std::streampos end = inputFile.tellg();

    	p->stream_size = (end - begin) / sizeof(unsigned long long int);
    	inputFile.seekg(0, std::ios::beg);

        stream = (unsigned long long int *) calloc(p->stream_size, sizeof(unsigned long long int));

    	inputFile.read((char*)stream, p->stream_size * sizeof(unsigned long long int));

    	inputFile.close();
    	
    	for (long i = 0; i < p->stream_size; i++){ 
    		unsigned int y = stream[i] >> 32;
			unsigned long long int x = stream[i] & 0x00000000FFFFFFFF;
			
			stream[i] = (x<<32) | y;
		}
    	
    	std::string paramsFileName = std::string(p->testFile) + "-params";

	    inputFile.open(paramsFileName, std::ios::in);
    	if (!inputFile.is_open()) {
        	std::cerr << "Error: bad parameters file" << std::endl;
        	exit (1);
    	}
    
    	std::string line;
    	getline(inputFile, line);
    	if (inputFile.eof())
        	exit (1);
    	getline(inputFile, line);
    	p->seed = (uint32_t)stol(line);
    	if (inputFile.eof())
        	exit (1);
	    getline(inputFile, line);
    	p->zipfpar = stod(line);
	    if (inputFile.eof())
        	exit (1);
	    getline(inputFile, line);
    	p->hwzpar = stod(line);
	    if (inputFile.eof())
        	exit (1);
	    getline(inputFile, line);
    	if (inputFile.eof())
        	exit (1);
	    getline(inputFile, line);
    	if (inputFile.eof())
        	exit (1);
	    getline(inputFile, line);
    	if (inputFile.eof())
        	exit (1);
    	getline(inputFile, line);
    	p->primary_distinct_items = p->secondary_distinct_items = (uint32_t)stol(line);
    	inputFile.close();

        free(p->testFile);
    }

    return(stream);
}

frequent_items* RunExact(unsigned long long int* s, input_parameters *p)
{
    int i;
   	unsigned int x, y;
	frequent_items *exact;

    exact = new frequent_items;
    
    unsigned int *xHist = (unsigned int*)calloc(p->primary_distinct_items+1, sizeof(unsigned int));
    unsigned int *yHist = (unsigned int*)calloc(p->secondary_distinct_items+1, sizeof(unsigned int));
    
    for (int i = 0; i < p->stream_size; i++) {
    	x = s[i] >> 32;
        xHist[x]++;
    }

    for (i = 1; i <= p->primary_distinct_items; i++) {
        if (xHist[i] > p->phi1 * p->stream_size) {
		    ChhCounter primarydata;
			
			primarydata.item = i;
			primarydata.count = xHist[i];
		    exact->primary[i] = primarydata;
		}
	}

    for (auto hhx : exact->primary) {
        for (int i = 0; i < p->stream_size; i++) {
        	x = s[i] >> 32;
	    	y = s[i] & 0x00000000FFFFFFFF;

            if (x == hhx.first)
                yHist[y]++;
        }

        for (int j = 1; j <= p->secondary_distinct_items; j++) {
            if (yHist[j] > p->phi2 * hhx.second.count) {
				
			    ChhCounter correlateddata;

				correlateddata.item = ((unsigned long long int)hhx.first) << 32 | j;
				correlateddata.count = yHist[j];	
				exact->correlated[correlateddata.item] = correlateddata;
            }
        }
        memset (yHist, 0, (p->secondary_distinct_items+1)* sizeof(unsigned int));
    }
    free(yHist);
    free(xHist);
                
    return exact;
}

unsigned int get_corr_exact_freq(input_parameters* p, unsigned long long int* s, unsigned long long int item)
{
	unsigned int count = 0;
	
	for (int i = 0; i < p->stream_size; i++) {
		if (s[i] == item)
			count ++;
	}
	return count;
}

unsigned int get_prim_exact_freq(input_parameters* p, unsigned long long int* s, unsigned int item)
{
	unsigned int count = 0;
	
	for (int i = 0; i < p->stream_size; i++) {
		unsigned int x = s[i] >> 32;
		if (x == item)
			count ++;
	}
	return count;
}


/* input parameters
 * input_parameters *p: configuration parameters
 * frequent_items *e: exact frequent correlated heavy hitters
 * frequent_items *chh: approximated frequent correlated heavy hitters
 * sketch_type *s: sketches data types (D1 and D*)
 * diagnostic *diag: diagnostic parameters
 */


void CheckOutput (input_parameters *p,unsigned long long int* s, frequent_items *e, frequent_items *chh, freq_type *freq, diagnostic *diag)
{
	std::map<unsigned long long int,ChhCounter>::iterator it;

	//Search true correlated heavy hitters
	
	diag->exact_num_primary = (int)e->primary.size();
   	diag->exact_num_correlated = (int)e->correlated.size();
	
	diag->chh_num_primary_candidates = (int)chh->primary.size();
	diag->chh_num_correlated_candidates = (int)chh->correlated.size();

	diag->chh_num_primary_exact = 0;
	diag->chh_num_correlated_exact = 0;
	
	//Computes the number of matches between the approximated chh and the exact chh 
    for (it = chh->primary.begin(); it != chh->primary.end(); ++it) {
		unsigned int x = (unsigned int)it->first;
		
		if (e->primary.find(x) != e->primary.end())
			diag->chh_num_primary_exact ++;
	}	
    for (it = chh->correlated.begin(); it != chh->correlated.end(); ++it) {
		unsigned long long int xy = it->first;
		
		if (e->correlated.find(xy) != e->correlated.end())
			diag->chh_num_correlated_exact ++;
	}	
	
	diag->recall_primary = (float)diag->chh_num_primary_exact/(float)diag->exact_num_primary;
	diag->precision_primary = (float)diag->chh_num_primary_exact/(float)diag->chh_num_primary_candidates;
	diag->recall_correlated = (float)diag->chh_num_correlated_exact/(float)diag->exact_num_correlated;
	diag->precision_correlated = (float)diag->chh_num_correlated_exact/(float)diag->chh_num_correlated_candidates;
	
	diag->ss_size = p->s1 * (sizeof(int) + sizeof(long));
	diag->ss_size += p->s1 * (p->s2 * (sizeof(int) + sizeof(long)) ) ;
	
	//Relative Error on Frequency estimation for Primary Items: Maximum, Minimum, Mean
	
	diag->max_rel_error_correlated_candidate = 0;
	diag->min_rel_error_correlated_candidate = 1;
	diag->mean_rel_error_correlated_candidate = 0;
	
	diag->max_abs_error_correlated_candidate = 0;
	diag->mean_abs_error_correlated_candidate = 0;

	diag->tolerance_eps1 = 0.0;
	diag->tolerance_eps2 = 0.0;

	//Evaluate the error for the frequency estimation

	for (auto chhx : chh->primary) {
		unsigned int x = (unsigned int)chhx.first;
		unsigned int exact_f_x;
		
		it = e->primary.find(x);
		if (it != e->primary.end())
			exact_f_x = it->second.count;
		else
		 	exact_f_x = get_prim_exact_freq(p, s, x);
							
		float tol_eps1 = p->phi1 - (float)exact_f_x / (float)p->stream_size;
		if (tol_eps1 > diag->tolerance_eps1)
			diag->tolerance_eps1 = tol_eps1;
	}


	for (auto chhxy : chh->correlated) {
		unsigned long long int xy = chhxy.first;
		int i = xy >> 32;
		int exact_f_x;
		
		int exact_f_xy;
		
		it = e->correlated.find(xy); 
		if (it != e->correlated.end())
			exact_f_xy = it->second.count;
		else
		 	exact_f_xy = get_corr_exact_freq(p, s, xy);

		int approximated_f_xy = chhxy.second.count;

		int abs_err_xy = exact_f_xy - approximated_f_xy;
		float rel_err_xy = (float)(abs_err_xy)/(float)exact_f_xy;
	
		if (rel_err_xy < 0) {
			fprintf(stderr,"Assertion failed: Frequency estimate %u for correlated (%d, %d) is less than the exact frequency %u (LINE: %d)\n", approximated_f_xy, i, (int)(xy & 0x00000000FFFFFFFF), exact_f_xy, __LINE__);
			exit (1);
		}
		
		if (rel_err_xy > diag->max_rel_error_correlated_candidate)
			diag->max_rel_error_correlated_candidate = rel_err_xy;
		if ( (approximated_f_xy != exact_f_xy) && rel_err_xy < diag->min_rel_error_correlated_candidate)
			diag->min_rel_error_correlated_candidate = rel_err_xy;			
		diag->mean_rel_error_correlated_candidate += rel_err_xy;

		if (abs_err_xy > diag->max_abs_error_correlated_candidate)
			diag->max_abs_error_correlated_candidate = abs_err_xy;
		diag->mean_abs_error_correlated_candidate += abs_err_xy;
		
		

		it = e->primary.find(i); 
		if (it != e->primary.end())
			exact_f_x = it->second.count;
		else
		 	exact_f_x = get_prim_exact_freq(p, s, i);
				
		float tol_eps2 = p->phi2 - (float)exact_f_xy / (float)exact_f_x;
		if (tol_eps2 > diag->tolerance_eps2)
			diag->tolerance_eps2 = tol_eps2;

	}	
	diag->mean_rel_error_correlated_candidate = diag->mean_rel_error_correlated_candidate / (float)diag->chh_num_correlated_candidates;
	diag->mean_abs_error_correlated_candidate = diag->mean_abs_error_correlated_candidate / (float)diag->chh_num_correlated_candidates;
		
}


void PrintDiagnostics(input_parameters *p, diagnostic *diag)
{
	std::string color_prim = "", color_corr = "";
	
	if (diag->recall_primary < 1)
		color_prim = "\x1b[31m";
	if (diag->recall_correlated < 1)
		color_corr = "\x1b[31m";
		
	
	printf("Frequent primary:  candidates (C) = %d \t true positives (P) = %d \t total exacts (E) = %d\n", diag->chh_num_primary_candidates, diag->chh_num_primary_exact, diag->exact_num_primary);
	printf("                   %s recall (P/E) = %f \x1b[0m \t precision (P/C) = %f\n",color_prim.c_str(), diag->recall_primary, diag->precision_primary);
			
	printf("Frequent correlated:  candidates (C) = %d \t true positives (P) = %d \t total exacts (E) = %d\n",diag->chh_num_correlated_candidates, diag->chh_num_correlated_exact, diag->exact_num_correlated);
	printf("                      %s recall (P/E) = %f \x1b[0m \t precision (P/C) = %f\n", color_corr.c_str(), diag->recall_correlated, diag->precision_correlated);
	printf("Relative Error on Correlated Candidate: Max = %f \t Min = %f \t Mean = %f \n", diag->max_rel_error_correlated_candidate, diag->min_rel_error_correlated_candidate, diag->mean_rel_error_correlated_candidate);
	printf("Absolute Error on Correlated Candidate: Max = %f \t Mean = %f \n\n", diag->max_abs_error_correlated_candidate, diag->mean_abs_error_correlated_candidate);
	
	printf("Measured values of: eps1 = %f (<= %f) \t eps2 = %f (<= %f) \n", diag->tolerance_eps1, p->eps1, diag->tolerance_eps2, p->eps2);
	
	printf("Execution time: Update_per_item = %f ms. \t Query = %f ms. \t Update_per_ms = %f\n", diag->update_time, diag->query_time, 1.0/diag->update_time);

	printf("Memory footprint: %.2f KB\n", (float)diag->ss_size/1024.0);
	
	//CSV output sent to the file descriptor 3

	//File descriptor 3 could have been opened by the shell due to I/O redirection	
	FILE* fp = fdopen(3, "w");
    if (!fp){
    	//File descriptor 3 has not been opened. We duplicate the stdout
    	dup2(1, 3);
    	fp = fdopen(3,"w");
    	printf("CSV Output:\n");
    }
    else
    	printf("CSV Output redirected to file\n");
    
    fprintf(stdout, "ni," \
					"Distinct_Prim," \
					"Distinct_Corr," \
					"seed," \
					"hseed," \
					"zipfpar," \
					"hwzpar," \
					"phi1," \
					"phi2," \
					"eps1," \
					"eps2," \
					"TrueHH_Prim," \
					"TrueHH_Corr," \
					"Cand_Prim," \
					"Cand_Corr," \
					"MatchHH_Prim," \
					"MatchHH_Corr," \
					"Recall_Prim," \
					"Precision_Prim," \
					"Recall_Corr," \
					"Precision_Corr," \
					"meanAbsErrCan_Corr," \
					"maxAbsErrCan_Corr," \
					"meanRelErrCan_Corr," \
					"maxRelErrCan_Corr," \
					"Measured_eps1," \
					"Measured_eps2," \
					"updatePerMs," \
					"queryTime(ms)," \
					"space(Kb)," \
					"s1," \
					"s2\n");

    fprintf(stderr, "%lu,", p->stream_size);
    fprintf(stderr, "%d,", p->primary_distinct_items);
    fprintf(stderr, "%d,", p->secondary_distinct_items);
    fprintf(stderr, "%d,", p->seed);
    fprintf(stderr, "%d,", p->hseed);
    fprintf(stderr, "%.2f,", p->zipfpar);
    fprintf(stderr, "%.2f,", p->hwzpar);
    fprintf(stderr, "%.3f,", p->phi1);
    fprintf(stderr, "%.3f,", p->phi2);
	fprintf(stderr, "%f,", p->eps1);
	fprintf(stderr, "%f,", p->eps2);
	fprintf(stderr, "%d,", diag->exact_num_primary);
	fprintf(stderr, "%d,", diag->exact_num_correlated);
	fprintf(stderr, "%d,", diag->chh_num_primary_candidates);
	fprintf(stderr, "%d,", diag->chh_num_correlated_candidates);
	fprintf(stderr, "%d,", diag->chh_num_primary_exact);
	fprintf(stderr, "%d,", diag->chh_num_correlated_exact);
	fprintf(stderr, "%f,", diag->recall_primary);
	fprintf(stderr, "%f,", diag->precision_primary);
	fprintf(stderr, "%f,", diag->recall_correlated);
	fprintf(stderr, "%f,", diag->precision_correlated);

	fprintf(stderr, "%f,", diag->mean_abs_error_correlated_candidate);
	fprintf(stderr, "%f,", diag->max_abs_error_correlated_candidate);
	fprintf(stderr, "%f,", diag->mean_rel_error_correlated_candidate);
	fprintf(stderr, "%f,", diag->max_rel_error_correlated_candidate);

	fprintf(stderr, "%f,", diag->tolerance_eps1);
	fprintf(stderr, "%f,", diag->tolerance_eps2);

	fprintf(stderr, "%f,", 1.0/diag->update_time);
	fprintf(stderr, "%f,", diag->query_time);

	fprintf(stderr, "%f,", (float)diag->ss_size/1024.0);

	fprintf(stderr, "%d,", p->s1);
	fprintf(stderr, "%d", p->s2);

	fprintf(stderr,"\n");
    fclose(fp);	
	
}



/*************************************************
 * Functions for handling frequent data structure
 *************************************************/


void frequent_Print(frequent_items *fi)
{
	map<unsigned long long int,ChhCounter>::iterator pf_it;
	printf("Number of Primary Frequent Items: %lu\n",  fi->primary.size());
    for (pf_it = fi->primary.begin(); pf_it != fi->primary.end(); ++pf_it) {
		ChhCounter sslc;
		sslc = pf_it->second;
		printf("<%u (%u)>, ",(unsigned int)sslc.item, sslc.count); 
	}
	printf("\n");

	printf("Number of Correlated Frequent Items: %lu\n",  fi->correlated.size());
    for (pf_it = fi->correlated.begin(); pf_it != fi->correlated.end(); ++pf_it) {
		unsigned int x, y;
		ChhCounter sslc;
		sslc = pf_it->second;
		x = sslc.item >> 32;
		y = sslc.item & 0x00000000FFFFFFFF;
		printf("<%u %u (%u)>, ",x, y, sslc.count); 
	}
	printf("\n");
}
