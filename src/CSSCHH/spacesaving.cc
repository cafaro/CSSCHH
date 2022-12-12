/**************************************************************************
 Implementation of SpaceSaving algorithm to Find Frequent Items
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
 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "spacesaving.h"

SSL_type * SSL_Init(long k)
{
    long i;

    SSL_type* result = (SSL_type*) calloc(1,sizeof(SSL_type));

    result->a= (long long) 698124007;
    result->b= (long long) 5125833;
    if (k<2) k=2;
    result->k=k;
    result->n=0;

    result->tblsz=SSL_HASHMULT*k;
    result->hashtable=(SSLCOUNTER**) calloc(result->tblsz,sizeof(SSLCOUNTER *));
    result->groups=(SSLGROUP *)calloc(k,sizeof(SSLGROUP));
    result->counters=(SSLCOUNTER *) calloc(k,sizeof(SSLCOUNTER));
    result->freegroups=(SSLGROUP **) calloc(k,sizeof(SSLGROUP*));

    for (i=0; i<result->tblsz; i++)
        result->hashtable[i]=NULL;

    result->root=result->groups;
    result->groups->count=0;
    result->groups->nextg=NULL;
    result->groups->previousg=NULL;

    result->groups->counters=result->counters;
    for (i=0; i<k; i++)
        result->freegroups[i]=&result->groups[i];
    result->gpt=1; // initialize list of free groups

    for (i=0; i<k; i++) {
        result->counters[i].item=0;
        result->counters[i].lerror=0;
        result->counters[i].hash=0;
        result->counters[i].nexti=NULL;
        result->counters[i].previousi=NULL;  // initialize values

        result->counters[i].parentg=result->groups;
        result->counters[i].nexting=&(result->counters[i+1]);
        result->counters[i].previousing=&(result->counters[i-1]);
        // create doubly linked list
    }
    result->counters[0].previousing=&(result->counters[k-1]);
    result->counters[k-1].nexting=&(result->counters[0]);
    // fix start and end of linked list

    return(result);
}

// void SSL_ShowGroups(SSL_type * ssl)
// {
// 	SSLGROUP *g;
// 	SSLCOUNTER *i,*first;
// 	long n;
//     SSLItemcount wt;

// 	g=ssl->root;
// 	wt=0;
// 	n=0;
// 	while (g!=NULL)
// 	{
// 		printf("Group %d :",g->count);
// 		first=g->counters;
// 		i=first;
// 		if (i!=NULL)
// 			do
// 			{
// 				printf("%d -> ",i->item);
// 				i=i->nexting;
// 				wt+=g->count;
// 				n++;
// 			}
// 		while (i!=first);
// 		else printf(" empty");
// 		printf(")");
// 		g=g->nextg;
// 		if ((g!=NULL) && (g->previousg->nextg!=g))
// 			printf("Badly linked");
// 		printf("\n");
// 	}
// 	printf("In total, %d counters, with a total count of %d\n",n,wt);
// }

void SSL_InsertIntoHashtable(SSL_type *ssl, SSLCOUNTER *newi, long i, SSLItem newitem)
{
    newi->nexti=ssl->hashtable[i];
    newi->item=newitem; // overwrite the old item
    newi->hash=i;
    newi->previousi=NULL;
    // insert item into the hashtable
    if (ssl->hashtable[i])
        ssl->hashtable[i]->previousi=newi;
    ssl->hashtable[i]=newi;
}

int SSL_cmp( const void * a, const void * b)
{
    SSLCOUNTER * x = (SSLCOUNTER*) a;
    SSLCOUNTER * y = (SSLCOUNTER*) b;
    if (x->parentg->count<y->parentg->count) return -1;
    else if (x->parentg->count>y->parentg->count) return 1;
    else return 0;
}

//std::map<uint32_t, uint32_t> SSL_Output(SSL_type * ssl, int thresh)
//{
//	std::map<uint32_t, uint32_t> res;
//
//	for (int i=0; i<ssl->k; ++i)
//		if (ssl->counters[i].parentg->count > thresh)
//			res.insert(std::pair<uint32_t, uint32_t>(ssl->counters[i].item, ssl->counters[i].parentg->count));
//
//	return res;
//}

long SSL_Output(SSL_type *ssl, SSLItemcount thresh, SSLCounterdata *counters)
{
    SSLGROUP *g;
    SSLCOUNTER *i,*first;
    long c;
    long z;

    g = ssl->root;
    c = 0;
    z = 0;
    if (thresh > 0) {
        while ((g != NULL) && (g->count < thresh)) {
            first = g->counters;
            i = first;
            if (i != NULL)
                do {
                    i=i->nexting;
                    z++;
                } while (i != first);
            g = g->nextg;
        }
    }

    while (g != NULL) {
        first = g->counters;
        i = first;
        if (i != NULL)
            do {
                counters[c+z].item = i->item;
                counters[c+z].lerror = i->lerror;
//                counters[c+z].rerror = i->rerror;
                counters[c+z].count = g->count;
                i=i->nexting;
                c++;
            } while (i != first);
        g = g->nextg;
    }
    return c;
}


SSLCOUNTER * SSL_GetNewCounter(SSL_type * ssl)
{
    SSLCOUNTER *newc;
    long j;

    newc = ssl->root->counters;  // take a counter from the first group
    // but currently it remains in the same group

    newc->nexting->previousing = newc->previousing;
    newc->previousing->nexting = newc->nexting;
    // unhook the new item from the linked list in the hash table

    // need to remove this item from the hashtable
    j = newc->hash;
    if (ssl->hashtable[j] == newc)
        ssl->hashtable[j] = newc->nexti;

    if (newc->nexti != NULL)
        newc->nexti->previousi = newc->previousi;
    if (newc->previousi != NULL)
        newc->previousi->nexti = newc->nexti;

    return (newc);
}

void SSL_PutInNewGroup(SSL_type * ssl, SSLCOUNTER *counter, SSLGROUP *newg)
{
    SSLGROUP * oldgroup;

    oldgroup = counter->parentg;
    // put item in the newg group
    counter->parentg = newg;

    if (counter->nexting != counter) { // if the group does not have size 1
        // remove the item from its current group
        counter->nexting->previousing = counter->previousing;
        counter->previousing->nexting = counter->nexting;
        oldgroup->counters = oldgroup->counters->nexting;
    } else { // group will be empty
        if (oldgroup->nextg != NULL) // there is another group
            oldgroup->nextg->previousg = oldgroup->previousg;
        if (ssl->root == oldgroup) // this is the first group
            ssl->root = oldgroup->nextg;
        else
            oldgroup->previousg->nextg=oldgroup->nextg;
        ssl->freegroups[--ssl->gpt]=oldgroup;
        // if we have created an empty group, remove it
    }
    counter->nexting=newg->counters;
    counter->previousing=newg->counters->previousing;
    counter->previousing->nexting=counter;
    counter->nexting->previousing=counter;
}

void SSL_AddNewGroupAfter(SSL_type * ssl, SSLCOUNTER *counter, SSLGROUP *oldgroup, SSLGROUP *after_this, SSLItemcount c)
{
    SSLGROUP *newgroup;

    // remove item from old group...
    counter->nexting->previousing = counter->previousing;
    counter->previousing->nexting = counter->nexting;
    oldgroup->counters = counter->nexting;

    //get new group
    newgroup = ssl->freegroups[ssl->gpt++];
    newgroup->count = oldgroup->count+c; // set count to be c more the prev group
    newgroup->counters = counter;

    newgroup->previousg = after_this;
    newgroup->nextg = after_this->nextg;
    if (after_this->nextg != NULL) // if there is another group
        after_this->nextg->previousg = newgroup;
    after_this->nextg = newgroup;

    counter->parentg = newgroup;
    counter->nexting = counter;
    counter->previousing = counter;
}

void SSL_MoveGroupAfter(SSL_type * ssl, SSLGROUP *to_move, SSLGROUP *after_this, SSLItemcount diffcount)
{
    to_move->count += diffcount;

    if (to_move != after_this) {
        // remove oldgroup from current position
        if (to_move->nextg != NULL)
            to_move->nextg->previousg = to_move->previousg;
        if (ssl->root == to_move)
            ssl->root = to_move->nextg;
        else
            to_move->previousg->nextg = to_move->nextg;

        // insert oldgroup in new position
        to_move->nextg = after_this->nextg;
        if (after_this->nextg != NULL)
            after_this->nextg->previousg = to_move;
        to_move->previousg = after_this;
        after_this->nextg = to_move;
    }
}


void SSL_IncrementCounter(SSL_type * ssl, SSLCOUNTER *counter, SSLItemcount c)
{
    SSLGROUP *oldgroup;
    SSLGROUP *tmpgrp;

    oldgroup = counter->parentg;
    tmpgrp = oldgroup;

    while ((tmpgrp->nextg != NULL) && (tmpgrp->nextg->count - oldgroup->count) < c) {
        tmpgrp = tmpgrp->nextg;
    }

    if ((tmpgrp->nextg != NULL) && (tmpgrp->nextg->count - oldgroup->count) == c) {
        SSL_PutInNewGroup(ssl,counter,tmpgrp->nextg);
    } else {
        if (counter->nexting == counter) // if there is only one item in the group...
            SSL_MoveGroupAfter(ssl,oldgroup,tmpgrp,c);
        else
            SSL_AddNewGroupAfter(ssl,counter,oldgroup,tmpgrp,c);
    }
}

void SSL_Update(SSL_type * ssl, SSLItem newitem, SSLItemcount c, SSLError lerror)
{
    if (c == 0) return; //nothing to do
    long h;
    SSLCOUNTER *il;

    ssl->n++;
    h=hash31(ssl->a,ssl->b,newitem) % ssl->tblsz;
    il=ssl->hashtable[h];
    while (il) {
        if (il->item ==newitem)
            break;
        il=il->nexti;
    }
    if (il==NULL) { // item is not monitored (not in hashtable)
        il=SSL_GetNewCounter(ssl);
//        /// and put it into the hashtable for the new item
        il->lerror = ssl->root->count + lerror;
        // initialize lerror with count of first group
        SSL_InsertIntoHashtable(ssl,il,h,newitem);
        // put the new counter into the first group
        // counter is already in first group by defn of how we got it
        SSL_IncrementCounter(ssl, il, c);
    } else {
        il->lerror += lerror;
        SSL_IncrementCounter(ssl, il, c);
        // if we have an item, we need to increment its counter
    }
}

unsigned int SSL_PointEstimate(SSL_type * ssl, SSLItem item)
{
    long h;
    SSLCOUNTER *il;

    ssl->n++;
    h=hash31(ssl->a,ssl->b,item) % ssl->tblsz;
    il=ssl->hashtable[h];
    while (il) {
        if (il->item == item)
            break;
        il=il->nexti;
    }
    if (il==NULL)
		return ssl->root->count;
	return il->parentg->count;
}

long SSL_Size(SSL_type * ssl)
{
    return sizeof(SSL_type)+(ssl->tblsz)*sizeof(SSLCOUNTER*) +
           (ssl->k)*(sizeof(SSLCOUNTER) + sizeof(SSLGROUP) + sizeof(SSLCOUNTER*));
}

void SSL_Destroy(SSL_type * ssl)
{
    free(ssl->freegroups);
    free(ssl->counters);
    free(ssl->groups);
    free(ssl->hashtable);
    free (ssl);
}
