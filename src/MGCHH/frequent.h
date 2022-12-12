//frequent.h -- simple frequent items routine
// see Misra&Gries 1982, Demaine et al 2002, Karp et al 2003
// implemented by Graham Cormode, 2002,2003

#include <vector>

struct FCounter {
    unsigned int item;
	void* value;
    int count;
};

typedef struct itemlist ITEMLIST;
typedef struct group GROUP;

struct group 
{
  int diff;
  ITEMLIST *items;
  GROUP *previousg, *nextg;
};

struct itemlist 
{
  int item;
  void* value;
  GROUP *parentg;
  ITEMLIST *previousi, *nexti;
  ITEMLIST *nexting, *previousing;
};

typedef struct freq_type{
  ITEMLIST **hashtable;
  GROUP *groups;
  ITEMLIST *items;
  int k;
  int tblsz;
  long long a,b;
  void (*free_value)(void*); //to be used to deallocate the value datastructure
} freq_type;

extern freq_type * Freq_Init(int num_counters, void (*free_value)(void*));
extern void Freq_Destroy(freq_type *freq);
extern ITEMLIST * Freq_Update(freq_type *, int newitem);
extern int Freq_Size(freq_type *freq);
//extern unsigned int *Freq_Output(freq_type *);
extern std::vector<FCounter> Freq_Output(freq_type *freq);
extern int Freq_PointEstimateCorr(freq_type* freq, int item_p, int item_s);
extern int Freq_PointEstimate(freq_type* freq, int item);
extern void SubtractCounter(ITEMLIST *newi);

