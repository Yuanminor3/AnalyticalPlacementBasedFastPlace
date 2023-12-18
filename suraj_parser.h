// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)

#ifndef __SURAJ_PARSER__H
#define __SURAJ_PARSER__H
#include<map>

using namespace std;

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

extern map <const char *, int, ltstr> nodeNameToNodeNum_map;

extern int numCellPins, // number of all terminals connected to the end points of (hyper) edges
	numhyper, 			// number of edges and hyperedges
	numCellsAndPads, 	// total number of movable cells (generall with names starting with a), and I/O pads (generally starting with p)
	numCells_noPads;	// total number of movable cells

extern int *cellPinArray;		// array of cell names used as endpoints of (hyper) edges
						// The size of the array is numCellPins

extern int *hEdge_idxToFirstEntryInPinArray;
						// (hyper)edge definitions.
						// The ith entry of this array is the
						// index (=j) of the first pin used in the (hyper)edge
						// The index j means the first pin of the (hyper) edge is the jth
						// entry in cellPinArray.
						// The degree of (hyper)edge i is determined
						// by subtracting the value of element i from the value of element i+1
						// The size of the array is numhyper+1

extern int *hyperwts;			// (hyper) edge weights. The size of the array is numhyper
extern int *vertexSize;		// cell and I/O pad sizes. The size of the array is numCellsAndPads

struct SPinLocation {
	int x, y;
};

extern SPinLocation *pinLocations;
						// The x,y location of the pins. 
						// The size of the array is (numCellsAndPads - numCells_noPads)
						// IMPORTANT: the 0th entry in this array, corresponds to p0, which has the cell index numCells_noPads
						// Similarly, pinLocations[i] corresponds to the (x,y) location of cell i+numCells_noPads.
						// Conversely, if you want to get the (x,y) location of cell j -- assuming j is the cell number associated with an 
						// I/O pad -- then j will be the p = (j - numCells_noPads)th I/O pad. pinLocations[p] has the (x,y) coordinate
						// you have to 


int parseIbmFile(char *inareFileName, char *innetFileName, char *inPadLocationFileName);


#endif 