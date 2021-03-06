#ifndef __VH_CREATE_CLUSTER_DELETION__H__
#define __VH_CREATE_CLUSTER_DELETION__H__

#include "vh_divethandler.h"
#include "vh_heap.h"
#include "../processbam.h"

#define LEFTSIDE 0
#define RIGHTSIDE 1

typedef struct MappingOnGenome
{
	DivetRow *readMappingPtr;
	struct MappingOnGenome *next;
} MappingOnGenome;


typedef struct RightBrkPointInterval
{
	int locBrkPointLeft;		// the left end of the right breakpoint interval
	int locBrkPointRight;		// the right end of the right breakpoint interval
	int key;			// Key is the locBrkPointLeft or locBrkPointRight
	char keyLorR;			// this field indicates that the key is the locBrkPointLeft or locBrkPointRight
	// the values it can get is 'L' or 'R'
	struct DivetRow *readMappingPtr;
} RightBrkPointInterval;

typedef struct outputElement{
	struct DivetRow *readMappingPtr;
	struct outputElement *next;
}outputElement;

typedef struct ClustersFound
{
	int *readMappingIdArray;
	struct DivetRow **readMappingPtrArray;
	int clusterSize;
	int leftBrkPoint;
	int isMaximalCluster;
	struct ClustersFound *next;
} ClustersFound;


void vh_initializeReadMapping_Deletion (sonic *this_sonic, int chr_index);
void vh_initializeReadMapping_Inversion (sonic *this_sonic, int chr_index);
void vh_initializeReadMapping_Insertion (sonic *this_sonic, int chr_index);
void vh_initializeReadMapping_TDup (sonic *this_sonic, int chr_index);
void vh_initializeReadMapping_InvDup( sonic *this_sonic, int chr_index, int invdup_location);
void vh_initializeReadMapping_InterDup( sonic *this_sonic, int chr_index, int invdup_location);
void vh_finalizeReadMapping (char *, int);
void vh_finalizeReadMapping_InvDup (char *chromosome_name, int chroSize);
void vh_finalizeReadMapping_InterDup (char *chromosome_name, int chroSize);
void vh_createDeletionClusters (int);
void vh_createInversionClusters (int);
void vh_createInsertionClusters (int);
void vh_createTDupClusters (int);
void vh_createInvDupClusters (int chroSize, int);
void vh_createInterDupClusters (int chroSize, int);
int vh_compare (const void *, const void *);
int vh_compare2 (const void *, const void *);
int vh_compareReadName (const void *, const void *);
int vh_compareInt (const void *, const void *);
int vh_comparePtr (const void *, const void *);
int vh_noGap (char *, int, int);
int vh_max (int, int);
int vh_min (int, int);
void vh_freeLinkedList (struct MappingOnGenome *);
int vh_freeLinkList2(PEAlistEls *ptr);
void vh_finalizeReadMapping (char *, int);
void vh_copyElBrkPointIntr (int, int);
int vh_isItSubset (int *, int, int *, int);
int vh_outputCluster (struct ClustersFound *, char);
void vh_addToPotentialOutput (int, struct Heap *, char);
void vh_flushOut (int, char);
void vh_createIntersectingIntervals (int, char);
int vh_notBothDirections (struct ClustersFound *);

extern struct MappingOnGenome **g_genomeIndexStart;
extern struct RightBrkPointInterval *g_listRightBrkPointIntr;
extern struct RightBrkPointInterval *g_tempListRightBrkPointIntr;
extern int g_maxListBrkPointIntr;
extern int g_listRightBrkPointIntrCount;
extern int g_maxDeltaAmongLibs;
extern int max_chromosome_size;
extern struct ClustersFound *g_listPotClusterFound;

#endif
