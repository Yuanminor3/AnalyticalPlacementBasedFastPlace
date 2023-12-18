// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)


# include<iostream>
# include<stdio.h>
# include<map>
# include<vector>
#include <string.h>
#include <fstream>
#include <string>


using namespace std;

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

map <const char *, int, ltstr> nodeNameToNodeNum_map;

int numCellPins, 		// number of all terminals connected to the end points of (hyper) edges
	numhyper, 			// number of edges and hyperedges
	numCellsAndPads, 	// total number of movable cells (generall with names starting with a), and I/O pads (generally starting with p)
	numCells_noPads;	// total number of movable cells

int *cellPinArray;		// array of cell names used as endpoints of (hyper) edges
						// The size of the array is numCellPins

int *hEdge_idxToFirstEntryInPinArray;
						// (hyper)edge definitions.
						// The ith entry of this array is the
						// index (=j) of the first pin used in the (hyper)edge
						// The index j means the first pin of the (hyper) edge is the jth
						// entry in cellPinArray.
						// The degree of (hyper)edge i is determined
						// by subtracting the value of element i from the value of element i+1
						// The size of the array is numhyper+1

int *hyperwts;			// (hyper) edge weights. The size of the array is numhyper
int *vertexSize;		// cell and I/O pad sizes. The size of the array is numCellsAndPads

struct SPinLocation {
	int x, y;
};

SPinLocation *pinLocations;
						// The x,y location of the pins. 
						// The size of the array is (numCellsAndPads - numCells_noPads)
						// IMPORTANT: the 0th entry in this array, corresponds to p0, which has the cell index numCells_noPads
						// Similarly, pinLocations[i] corresponds to the (x,y) location of cell i+numCells_noPads.
						// Conversely, if you want to get the (x,y) location of cell j -- assuming j is the cell number associated with an 
						// I/O pad -- then j will be the p = (j - numCells_noPads)th I/O pad. pinLocations[p] has the (x,y) coordinate
						// you have to 

int parseIbmFile(char *inareFileName, char *innetFileName, char *inPadLocationFileName)
{
	char line[80];
	char vertexname[10],nodedesc[10];
	char isStartHyperNet;
	int hyperweight;
	FILE *innetFile, *inareFile, *inPadLocationFile;

	innetFile = fopen (innetFileName, "r");
        if (!innetFile) {
                printf ("ERROR: Cannot open input file %s.\n", innetFileName);
                return -1;
        }
	inareFile = fopen (inareFileName, "r");
        if (!inareFile) {
                printf ("ERROR: Cannot open input file %s.\n", inareFileName);
                return -1;
        }
	inPadLocationFile  = fopen (inPadLocationFileName, "r");
        if (!inPadLocationFile) {
                printf ("ERROR: Cannot open input file %s.\n", inPadLocationFileName);
                return -1;
        }
	
	

	fscanf(innetFile,"%*d %d %d %d %d\n",&numCellPins,&numhyper,&numCellsAndPads,&numCells_noPads);
	++numCells_noPads;

	cout << "numCellPins, numhyper, numCellsAndPads, numCells_noPads = " << numCellPins << ", " << numhyper << ", " << numCellsAndPads << ", " << numCells_noPads << endl;
	// numCellPins is the total number of end-points (pins) for hyperedges and edges.
	// numhyper is the total number of hyperedges + edges
	// the other two variables should be self-explanatory

	// For example, if we have this netlist:
	// Four vertices a0, a1, a2, p1
	// two edges (a0, a1), (a1, a2)
	// one hyperedge (p1, a0, a2)
	// Then numCellsAndPads == 4
	// numCells_noPads == 3
	// numhyper == 3
	// numCellPins = 7   (corresponding to how many times cell names appear
	//					  in the description of (hyper)edges )


	// Read cell names, create cell name --> cell number map
	vertexSize = (int *)malloc(numCellsAndPads*sizeof(int));

	for(int i=0; i<numCellsAndPads; i++)
	{
		fscanf(inareFile,"%s %d",vertexname,&vertexSize[i]);
		//printf("\n%s %d",vertexname,vertexSize[i]);
		nodeNameToNodeNum_map[strdup(vertexname)]=i;
		// if cell name starts with "a", it's a movable cell
		// and if starts with "p", it is an I/O pad.
	}


	// Read (hyper)edges
	hEdge_idxToFirstEntryInPinArray = (int *) malloc((numhyper+1)*sizeof(int));
	cellPinArray = (int *) malloc(numCellPins*sizeof(int));
	hyperwts = (int *) malloc(numhyper*sizeof(int));
	if(hEdge_idxToFirstEntryInPinArray==NULL || cellPinArray == NULL || hyperwts == NULL)
	{
		printf("\nUnable to allocate memory");
		return -1;
	}
	int hypercount = 0;
	int pinCount = 0;
	
	while (fgets(line, 80, innetFile) != NULL) 
   	{
		sscanf(line,"%s %c %d",nodedesc,&isStartHyperNet,&hyperweight);
		if(isStartHyperNet=='s')
		{
			hEdge_idxToFirstEntryInPinArray[hypercount] = pinCount;
			hyperwts[hypercount] = hyperweight;
			++hypercount;
		}
		else
		{
			sscanf(line,"%s %d",nodedesc,&hyperweight);
		}
		cellPinArray[pinCount] = nodeNameToNodeNum_map[nodedesc]; 
		++pinCount;
	}
	hEdge_idxToFirstEntryInPinArray[hypercount] = numCellPins;

	// cellPinArray lists cell indices involved in (hyper)edges
	// assuming hEdge_idxToFirstEntryInPinArray[i] = j, it means
	// that (hyper)edge i's first connected cell is cellPinArray[j]
	// its second cell is cellPinArray[j+1], and if it's a hyperedge,
	// its third cell is cellPinArray[j+2], and so on.
	// We know the degree K of (hyper)edge i by looking at
	// K = hEdge_idxToFirstEntryInPinArray[i+1] - hEdge_idxToFirstEntryInPinArray[i]

	// In the small graph example above, 
	// cellPinArray[] = {// this index 0 is the beginning of the first (hyper)edge #0
	//					 0, // corresponds to a0 in (a0, a1)
	//					 1, // corresponds to a1 in (a0, a1)
	//					// this index 2 is the beginning of (hyper) edge #1
	//					1,  // corresponds to a1 in the second edge (a1, a2)
	//					2,  // corresponds to a2 in the second edge (a1, a2)
	// 					// this index 4 is the beginning of the third hyperedge  (p1, a0, a2)
	//					3,  // corresponds to p1 in the third (hyper) edge (p1, a0, a2)
	//					0,  // corresponds to a1 in the third (hyper) edge (p1, a0, a2)
	//					2  // corresponds to a1 in the third (hyper) edge (p1, a0, a2)
	//					}

	// hEdge_idxToFirstEntryInPinArray[] = {0, 2, 4, 7}
	// The very last element is numCellPins, so that you can get
	// the degree of the last hyperedge by subtracting 4 from 7
	// 7-4 = 3.



	// Read I/O Pad locations
	int numPads = numCellsAndPads - numCells_noPads;
	pinLocations = new SPinLocation[numPads];
	for(int i=0; i<numPads; i++)
		{
			char padName[80];
			int x, y;
			fscanf(inPadLocationFile,"%s %d %d",padName,&x, &y);
			int padCellIdx = nodeNameToNodeNum_map[padName];
			pinLocations[padCellIdx-numCells_noPads].x = x;		// pay attention how it's indexed 
			pinLocations[padCellIdx-numCells_noPads].y = y;
		}


	fclose(inPadLocationFile);
	fclose(innetFile);
    fclose(inareFile);

	return 0;
}
		
