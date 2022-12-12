/**************************************************************************
 Implementation of SpaceSaving algorithm to Find Frequent Items in parallel
 Based on a paper by:
 Metwally, Agrawal, El Abbadi 2005
 Original sequential code by G. Cormode 2002, 2003

 Original Code: 2002-11, 2003-10
 This version: 2013-05

 This work is licensed under the Creative Commons
 Attribution-NonCommercial License. To view a copy of this license,
 visit http://creativecommons.org/licenses/by-nc/1.0/ or send a letter
 to Creative Commons, 559 Nathan Abbott Way, Stanford, California
 94305, USA.
 **************************************************************************/

#ifndef SPACESAVING_h
#define SPACESAVING_h

#include "prng.h"

#define SSL_HASHMULT 3

typedef uint64_t SSLItem;
typedef long     SSLItemcount;
typedef long     SSLError;

typedef struct ssl_counter SSLCOUNTER;
typedef struct ssl_group   SSLGROUP;

struct ssl_group {
    SSLItemcount count;
    SSLCOUNTER	 *counters;
    SSLGROUP	 *previousg, *nextg;
};

struct ssl_counter {
    SSLItem      item;
    SSLError     lerror;
    long	 hash;
    SSLGROUP     *parentg;
    SSLCOUNTER   *previousi, *nexti;
    SSLCOUNTER   *nexting, *previousing;
};

typedef struct SSL_type {
    SSLItemcount n;
    long         gpt;
    long         k;
    long         tblsz;
    long long    a,b;
    SSLGROUP     *root;
    SSLCOUNTER   *counters;
    SSLGROUP     *groups;
    SSLGROUP     **freegroups;
    SSLCOUNTER   **hashtable;
} SSL_type;

typedef struct SSLCounterdata {
    SSLItem      item;
    SSLItemcount count;
    SSLError     lerror;
    SSLError     rerror;
} SSLCounterdata;

extern SSL_type *SSL_Init(long);
extern void SSL_Update(SSL_type *, SSLItem, SSLItemcount, SSLError lerror);
extern long SSL_Output(SSL_type *,SSLItemcount, SSLCounterdata *);
extern void SSL_Destroy(SSL_type *);
//extern void SSL_ShowGroups(SSL_type *);
extern long SSL_Size(SSL_type *);
unsigned int SSL_PointEstimate(SSL_type * ssl, SSLItem item);


#endif


