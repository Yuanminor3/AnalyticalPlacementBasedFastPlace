// written by Suraj Dalvi, ECE Dept., Univ. of Minnesota
// Heavily modified by Kia Bazargan (renaming variables, adding
// comments, creating a clean version of the hyperedge data structure)

#include <fstream>
#include <string>
#include<iostream>
#include <algorithm> // For std::max_element and std::lower_bound
#include<stdio.h>
#include<map>
#include <utility>
#include <cmath>
#include<vector>
#include <string.h>
#include "suraj_parser.h"
#include <iomanip> // For formatting output
#include <sstream> // For stringstream



using namespace std;
double** QMatrix = nullptr; // Global pointer for the Q matrix
int qMatrixSize;
vector<double> xd; // Initialize xd with size qMatrixSize and all values set to 0
vector<double> yd; // Initialize yd with size qMatrixSize and all values set to 0
vector<double> solutionX; // Save X coordinates by Conjugate Gradient
vector<double> solutionY; // Save Y coordinates by Conjugate Gradient
vector<double> newXPositions;
vector<double> newYPositions;
double maxX, maxY;
double delta = 1.5;
int numBinsX, numBinsY;
vector<double> OBx, OBy;  // Boundaries of bins in regular structure
vector<double> NBx, NBy;  // Boundaries of bins in new, uneven structure
vector<int> Ux, Uy;       // Utilization of bins
vector< vector<int> > cellsThatBelongToBinX, cellsThatBelongToBinY;

float preWL; // Before spreading WL
float afterWL; // After spreading WL


// Function to print map (debug)
void printMap(const map<const char*, int, ltstr>& nodeNameToNodeNum_map) {
    for (const auto& pair : nodeNameToNodeNum_map) {
        cout << "Key: " << pair.first << " Value: " << pair.second << endl;
    }
}
// Function to print hEdge_idxToFirstEntryInPinArray (debug)
void printHEdgeIdxArray(int* array, int size) {
    for (int i = 0; i < size; ++i) {
        cout << "hEdge_idxToFirstEntryInPinArray[" << i << "] = " << array[i] << endl;
    }
}
// Function to print cellPinArray (debug)
void printCellPinArray(int* array, int size) {
    for (int i = 0; i < size; ++i) {
        cout << "cellPinArray[" << i << "] = " << array[i] << endl;
    }
}
// Function to print hyperwts array (debug)
void printHyperwts(int* array, int size) {
    for (int i = 0; i < size; ++i) {
        cout << "hyperwts[" << i << "] = " << array[i] << endl;
    }
}
// Function to print vertexSize array (debug)
void printVertexSize(int* array, int size) {
    for (int i = 0; i < size; ++i) {
        cout << "vertexSize[" << i << "] = " << array[i] << endl;
    }
}
// Function to print vertexSize array (debug)
void printPinLocations() {
     // Calculate the size of the pinLocations array
    int size = numCellsAndPads - numCells_noPads;
    // Print the (x, y) location of each pin
    for (int i = 0; i < size; ++i) {
        cout << "Pin " << i+1 << ": (" 
                  << pinLocations[i].x << ", " << pinLocations[i].y << ")" << endl;
    }
}
// Function to print xd vector and yd vector (debug)
void printXdYd() {
    stringstream ss;
    ss << "xd vector: [";
    for (size_t i = 0; i < qMatrixSize; ++i) {
        ss << xd[i];
        if (i < qMatrixSize - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    cout << ss.str() << endl;

	ss.str("");
	ss.clear();

    ss << "yd vector: [";
    for (size_t i = 0; i < qMatrixSize; ++i) {
        ss << yd[i];
        if (i < qMatrixSize - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    cout << ss.str() << endl;
}
// Function to print xd vector and yd vector (debug)
void printXYCor() {
    stringstream ss;
    ss << "x cordinates: [";
    for (size_t i = 0; i < qMatrixSize; ++i) {
        ss << solutionX[i];
        if (i < qMatrixSize - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    cout << ss.str() << endl;

	ss.str("");
	ss.clear();

    ss << "y cordinates: [";
    for (size_t i = 0; i < qMatrixSize; ++i) {
        ss << solutionY[i];
        if (i < qMatrixSize - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    cout << ss.str() << endl;
}
void printXYCor_after(){
    stringstream ss;
    ss << "x cordinates: [";
    for (size_t i = 0; i < newXPositions.size(); ++i) {
        ss << newXPositions[i];
        if (i < newXPositions.size() - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    cout << ss.str() << endl;

	ss.str("");
	ss.clear();

    ss << "y cordinates: [";
    for (size_t i = 0; i < newYPositions.size(); ++i) {
        ss << newYPositions[i];
        if (i < newYPositions.size() - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    cout << ss.str() << endl;

}

// Function to determine the size of the Q matrix
int determineQMatrixSize(const int* hEdge_idxToFirstEntryInPinArray, int size, int numCells_noPads) {
    int countGreaterThanThree = 0;

    for (int i = 0; i < size - 1; ++i) {
        if (hEdge_idxToFirstEntryInPinArray[i + 1] - hEdge_idxToFirstEntryInPinArray[i] > 3) {
            countGreaterThanThree++;
        }
    }

    return numCells_noPads + countGreaterThanThree;
}
// Function to print the Q matrix
void printQMatrix() {
	cout << "\nQ matrix:" << endl;
    for (int i = 0; i < qMatrixSize; ++i) {
        for (int j = 0; j < qMatrixSize; ++j) {
            cout << setw(8) << QMatrix[i][j] << " ";
        }
        cout << endl; // New line after each row
    }
}
// Function to fill the Q matrix
void fillQMatrix() {
    for (int a = 0; a < numCells_noPads; ++a) { // Iterate over each line (row)
        for (int i = 0; i < numhyper; ++i) { // Iterate part by part
            // Determine the start and end of the current part
            int start = hEdge_idxToFirstEntryInPinArray[i];
            int end = hEdge_idxToFirstEntryInPinArray[i + 1];

            bool containsLineA = false;
            int k = 0; // Number of elements in the part

            // Check if the part contains line number a and count elements
            for (int j = start; j < end; ++j) {
                if (cellPinArray[j] == a) {
                    containsLineA = true;
                }
                k++;
            }

            if (!containsLineA) continue; // Skip if the part doesn't contain line number a

            // Calculate the value to add to Q[a][a]
            if (k > 3) {
                QMatrix[a][a] += k / (k - 1.0);
            } else {
                double valueToAdd = 1.0 / (k - 1.0);
				if (k==3){
					QMatrix[a][a] += 1/(k - 1.0)*2;
				}else{
					QMatrix[a][a] += valueToAdd;

				}

                // Iterate each element in this part
                for (int j = start; j < end; ++j) {
                    int b = cellPinArray[j];
                    if (b < numCells_noPads && b != a && b > a) {
                        QMatrix[a][b] -= valueToAdd;
                        QMatrix[b][a] -= valueToAdd; // Assuming the matrix is symmetric
                    }
                }
            }
        }
    }

	// Deal with dummy points in Q matrix
	int aCounter = numCells_noPads; // Start from numCells_noPads for lines with more than 3 elements
    for (int i = 0; i < numhyper; ++i) { // Iterate part by part
        // Determine the start and end of the current part
        int start = hEdge_idxToFirstEntryInPinArray[i];
        int end = hEdge_idxToFirstEntryInPinArray[i + 1];

        int k = end - start; // Number of elements in the part
		
        if (k > 3) {
            int a = aCounter; // Corresponding line in Q matrix
            double value = k / static_cast<double>(k - 1);
			QMatrix[a][a] += k*k / static_cast<double>(k - 1);

            //Iterate all elements in the part
            for (int j = start; j < end; ++j) {
                int b = cellPinArray[j];
				//cout << "\n" << a << endl;

                if (b < numCells_noPads) {
                    QMatrix[a][b] -= value;
                    QMatrix[b][a] -= value;
                }
                
            }
            aCounter++; // Increment for the next line corresponding to parts with more than 3 elements
        }
    }
}
// Function to get xd and yd vector
void fillXdYd() {
	yd.resize(qMatrixSize, 0.0);
	xd.resize(qMatrixSize, 0.0);
	int s = numCells_noPads;
    for (int part = 0; part < numhyper; ++part) {
        int start = hEdge_idxToFirstEntryInPinArray[part];
        int end = hEdge_idxToFirstEntryInPinArray[part + 1];
		int k = end - start;
		if (k >= 4){
			for (int i = start; i < end; ++i) {
				int a = cellPinArray[i];
				if (a >= numCells_noPads) {
					double value = k / static_cast<double>(k - 1);
					xd[s] += pinLocations[a-numCells_noPads].x*value;
					yd[s] += pinLocations[a-numCells_noPads].y*value;
				}
			}
			
			s += 1;
			continue;
		}

        for (int i = start; i < end; ++i) {
            int a = cellPinArray[i];
            if (a >= numCells_noPads) {
                for (int j = start; j < end; ++j) {
                    int b = cellPinArray[j];
                    if (b < numCells_noPads) {
						if (k < 3){
							xd[b] += pinLocations[a-numCells_noPads].x;
							yd[b] += pinLocations[a-numCells_noPads].y;
						}else if(k==3){
							xd[b] += pinLocations[a-numCells_noPads].x*(0.5);
							yd[b] += pinLocations[a-numCells_noPads].y*(0.5);
						}else{
							continue;
						}
						
                    }
                }
            }
        }
    }
}
// Output to coordinates.txt file
void printCoordinatesToFile() {

    // Open file stream
    ofstream file("output/coordinates.txt");

    // Check if file is open
    if (!file.is_open()) {
        cerr << "Error opening file 'coordinates.txt'" << endl;
        return;
    }
	// matrix
	file << "\nQ matrix " << "(Size of " << qMatrixSize << "):" << endl;
    for (int i = 0; i < qMatrixSize; ++i) {
        for (int j = 0; j < qMatrixSize; ++j) {
            file << setw(8) << QMatrix[i][j] << " ";
        }
        file << endl; // New line after each row
    }

	file << endl;
	// xd and yd vectors
	stringstream ss;
    ss << "xd vector: [";
    for (size_t i = 0; i < qMatrixSize; ++i) {
        ss << xd[i];
        if (i < qMatrixSize - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    file << ss.str() << endl;

	ss.str("");
	ss.clear();

    ss << "yd vector: [";
    for (size_t i = 0; i < qMatrixSize; ++i) {
        ss << yd[i];
        if (i < qMatrixSize - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    file << ss.str() << endl;

	file << endl;
	file << "[**Before Spreading Coordinates**] Moveable Cells, Dummy node (stars), I/O Pads:\n" << endl;
    // Write coordinates
    for (int i = 0; i < qMatrixSize; ++i) {
        if (i < numCells_noPads) {
            file << "a" << i << " " << solutionX[i] << " " << solutionY[i] << "\n";
        } else{
            file << "d" << (i - numCells_noPads + 1) << " " << solutionX[i] << " " << solutionY[i] << "\n";
        }
    }

	// I/O Pads
	int size = numCellsAndPads - numCells_noPads;
    // Print the (x, y) location of each pin
    for (int i = 0; i < size; ++i) {
        file << "p" << i+1 << " " 
                  << pinLocations[i].x << " " << pinLocations[i].y << endl;
    }

	file << "\nSqrt of sum of square WL (only b/w movable cells) before spreading: " << preWL << endl;

    // Close file stream
    file.close();
}
// Output to coordinates2.txt file
void printCoordinatesToFile2() {

    // Open file stream
    ofstream file("output/coordinates2.txt");

    // Check if file is open
    if (!file.is_open()) {
        cerr << "Error opening file 'coordinates2.txt'" << endl;
        return;
    }
    
	file << "[**After Spreading Coordinates**] Moveable Cells, Dummy node (stars), I/O Pads:\n" << endl;
    // Write coordinates
    for (int i = 0; i < qMatrixSize; ++i) {
        if (i < numCells_noPads) {
            file << "a" << i << " " << newXPositions[i] << " " << newYPositions[i] << "\n";
        } else{
            file << "d" << (i - numCells_noPads + 1) << " " << newXPositions[i] << " " << newYPositions[i] << "\n";
        }
    }

	// I/O Pads
	int size = numCellsAndPads - numCells_noPads;
    // Print the (x, y) location of each pin
    for (int i = 0; i < size; ++i) {
        file << "p" << i+1 << " " 
                  << pinLocations[i].x << " " << pinLocations[i].y << endl;
    }

	file << "\nSqrt of sum of square WL (only b/w movable cells) after spreading: " << afterWL << endl;

    // Close file stream
    file.close();
}
// Conjugate Gradient
vector<double> conjugateGradient(const vector<double>& rhs) {
     vector<double> x(qMatrixSize, 0); // Solution vector initialized to 0
    vector<double> r = rhs; // Residual vector initialized to rhs
    vector<double> p = r; // The direction vector initialized to r
    vector<double> new_r(qMatrixSize);
    double alpha, beta, rnorm;

    for (int i = 0; i < qMatrixSize; ++i) {
        // Matrix-vector multiplication: QMatrix * p
        vector<double> Qp(qMatrixSize, 0);
        for (int j = 0; j < qMatrixSize; ++j) {
            for (int k = 0; k < qMatrixSize; ++k) {
                Qp[j] += QMatrix[j][k] * p[k];
            }
        }

        // Calculating alpha
        double dot_pQp = 0, dot_rr = 0;
        for (int j = 0; j < qMatrixSize; ++j) {
            dot_pQp += p[j] * Qp[j];
            dot_rr += r[j] * r[j];
        }
        alpha = dot_rr / dot_pQp;

        // Update x and r
        for (int j = 0; j < qMatrixSize; ++j) {
            x[j] += alpha * p[j];
            new_r[j] = r[j] - alpha * Qp[j];
        }

        // Check for convergence
        rnorm = 0;
        for (double val : new_r) {
            rnorm += val * val;
        }
        if (sqrt(rnorm) < 1e-10) {
            break;
        }

        // Calculating beta
        beta = 0;
        for (int j = 0; j < qMatrixSize; ++j) {
            beta += new_r[j] * new_r[j];
        }
        beta /= dot_rr;

        // Update p
        for (int j = 0; j < qMatrixSize; ++j) {
            p[j] = new_r[j] + beta * p[j];
        }

        r = new_r;
    }

    return x;
}
//Function to calculate WL before spreading
void calculateTotalWidthLength() {
    float totalDistance = 0.0f;
    float totalDistance2 = 0.0f;
    int dum = 0;
    for (int part = 0; part < numhyper; ++part) {
        int start = hEdge_idxToFirstEntryInPinArray[part];
        int end = hEdge_idxToFirstEntryInPinArray[part + 1];
        int k = end - start;
        int n = 0;
        vector<int> cels;
        for (int i = start; i < end; ++i) {
            if (cellPinArray[i] < numCells_noPads){
                cels.push_back(cellPinArray[i]);
                n++;
            }
        }        
        if (k>3){
            double w = k / static_cast<double>(k - 1);
            for (int i = 0; i < n; ++i){
                totalDistance += w*(pow(solutionX[cels[i]] - solutionX[numCells_noPads+dum], 2) + pow(solutionY[cels[i]] - solutionY[numCells_noPads+dum], 2));
                totalDistance2 += w*(pow(newXPositions[cels[i]] - newXPositions[numCells_noPads+dum], 2) + pow(newYPositions[cels[i]] - newYPositions[numCells_noPads+dum], 2));
            }
            dum += 1;

        }else if(k==3){
            double w = 0.5;
            for (int i = 0; i < n-1; ++i){
                for (int j = i+1; j < n; ++j){
                    totalDistance += w*(pow(solutionX[cels[i]] - solutionX[cels[j]], 2) + pow(solutionY[cels[i]] - solutionY[cels[j]], 2));
                    totalDistance2 += w*(pow(newXPositions[cels[i]] - newXPositions[cels[j]], 2) + pow(newYPositions[cels[i]] - newYPositions[cels[j]], 2));

                }
            }
        }else{
            double w=1.0;
            for (int i = 0; i < n-1; ++i){
                for (int j = i+1; j < n; ++j){
                    totalDistance += w*(pow(solutionX[cels[i]] - solutionX[cels[j]], 2) + pow(solutionY[cels[i]] - solutionY[cels[j]], 2));
                    totalDistance2 += w*(pow(newXPositions[cels[i]] - newXPositions[cels[j]], 2) + pow(newYPositions[cels[i]] - newYPositions[cels[j]], 2));
                }
            }
        }

    }
    preWL = sqrt(totalDistance);
    afterWL = sqrt(totalDistance2);

}
// Function to find the bin index for a given coordinate
int findBinIndex(const vector<double>& binBoundaries, double coordinate) {
    vector<double>::const_iterator it = lower_bound(binBoundaries.begin(), binBoundaries.end(), coordinate);
    int index = distance(binBoundaries.begin(), it);
    // Special case for when coordinate is exactly on a boundary
    if (index == 0) {
        index++;
    }
    return index;
}
// Assign cells to bins and calculate utilization
void assignCellsToBinsAndCalculateUtilization() {
    // Clear previous utilization and bin assignments
    Ux.assign(OBx.size(), 0);
    Uy.assign(OBy.size(), 0);
    cellsThatBelongToBinX.assign(OBx.size(), vector<int>());
    cellsThatBelongToBinY.assign(OBy.size(), vector<int>());

    // Assign cells to bins based on their initial positions and calculate utilization
    for (int i = 0; i < qMatrixSize; ++i) {
        int binIndexX = findBinIndex(OBx, solutionX[i]);
        int binIndexY = findBinIndex(OBy, solutionY[i]);

        cellsThatBelongToBinX[binIndexX].push_back(i);
        cellsThatBelongToBinY[binIndexY].push_back(i);
        Ux[binIndexX]++; // Since each cell is 1x1, increment the utilization by 1
        Uy[binIndexY]++;
    }
}
// Initializes the bin boundaries and utilizations
void initializeBins() {
    // Find the extents of the pin locations
    maxX = pinLocations[0].x;
    maxY = pinLocations[0].y;
    int numPins = numCellsAndPads - numCells_noPads;
    for (int i = 1; i < numPins; ++i) {
        if(maxX < pinLocations[i].x){
            maxX = pinLocations[i].x;
        }
        if(maxY < pinLocations[i].y){
            maxY = pinLocations[i].y;
        }
    }

    // Determine the number of bins on each side
    int numBinsSide = trunc(maxX/2.0);

    numBinsX = numBinsSide;
    numBinsY = numBinsSide; 

    // Calculate bin width and height
    double binWidth = maxX / numBinsX;
    double binHeight = maxY / numBinsY;

    // Initialize regular bin boundaries
    OBx.resize(numBinsX + 1);
    OBy.resize(numBinsY + 1);
    for (int i = 0; i <= numBinsX; ++i) {
        OBx[i] = i * binWidth;
    }
    for (int i = 0; i <= numBinsY; ++i) {
        OBy[i] = i * binHeight;

    }
    assignCellsToBinsAndCalculateUtilization();


    // Calculate NBx and NBy based on the utilizations and delta
    NBx.resize(OBx.size(),0);
    NBy.resize(OBy.size(),0);

    // Calculate NBx for the internal bins (excluding the first and last boundary)
    for (int i = 1; i < OBx.size() - 1; ++i) {
        NBx[i] = (OBx[i-1] * (Ux[i+1]/binWidth + delta) + OBx[i+1] * (Ux[i]/binWidth + delta)) / (Ux[i]/binWidth + Ux[i+1]/binWidth + 2 * delta);
    }

    // Repeat for NBy using OBy and Uy
    for (int i = 1; i < OBy.size() - 1; ++i) {
        NBy[i] = (OBy[i-1] * (Uy[i+1]/binHeight + delta) + OBy[i+1] * (Uy[i]/binHeight + delta)) / (Uy[i]/binHeight + Uy[i+1]/binHeight + 2 * delta);
    }

    // Set the boundaries for the first and last bins explicitly
    NBx[0] = OBx[0];
    NBx[OBx.size() - 1] = OBx[OBx.size() - 1];
    NBy[0] = OBy[0];
    NBy[OBy.size() - 1] = OBy[OBy.size() - 1];
}
// Function to calculate new cell positions based on shifting
pair< vector<double>, vector<double> > shiftCells() {
    // Initialize the new solution vectors with original positions
    vector<double> newSolutionX(solutionX);
    vector<double> newSolutionY(solutionY);

    // Iterate over each cell to compute its new position based on its bin
    for (int i = 0; i < qMatrixSize; ++i) {
        int binIndexX = findBinIndex(OBx, solutionX[i]); // Find the current bin index for cell i in the X direction
        int binIndexY = findBinIndex(OBy, solutionY[i]); // Find the current bin index for cell i in the Y direction
        double alpha_x = 0.8;
        double alpha_y = 0.8;


        newSolutionX[i] = (NBx[binIndexX]*(solutionX[i]-OBx[binIndexX-1])+NBx[binIndexX-1]*(OBx[binIndexX]-solutionX[i]))/(OBx[binIndexX]-OBx[binIndexX-1]);
        newSolutionX[i] = solutionX[i] + alpha_x * (newSolutionX[i] - solutionX[binIndexX]); // Apply alpha_x


        newSolutionY[i] = (NBy[binIndexY]*(solutionY[i]-OBy[binIndexY-1])+NBy[binIndexY-1]*(OBy[binIndexY]-solutionY[i]))/(OBy[binIndexY]-OBy[binIndexY-1]);
        newSolutionY[i] = solutionY[i] + alpha_y * (newSolutionY[i] - solutionY[binIndexY]); // Apply alpha_x
        
    }
    return make_pair(newSolutionX, newSolutionY);
}


// ***** Main Function *****	
int main(int argv, char *argc[])
{
	char inareFileName[100];
	char innetFileName[100];
	char inPadLocationFileName[100];

	if (argv!=2) {
		cout << "Please provide a circuit file name with no extension." << endl;
		return 1;
	}
        	     
    cout << "Reading circuit file " << argc[1] << endl;

	strcpy (inareFileName, argc[1]);
	strcat(inareFileName, ".are");
	strcpy(innetFileName,argc[1]);
	strcat(innetFileName,".net");
	strcpy(inPadLocationFileName,argc[1]);
	strcat(inPadLocationFileName,".kiaPad");

	int success = parseIbmFile(inareFileName, innetFileName, inPadLocationFileName);
	if (success == -1) {
		cout << "Error reading input file(s)" << endl;
		return 0;
	}

	printf("\nNumber of vertices,hyper = %d %d\n",numCellsAndPads,numhyper);

	// printMap(nodeNameToNodeNum_map);
	// printHEdgeIdxArray(hEdge_idxToFirstEntryInPinArray, numhyper + 1);
    // printCellPinArray(cellPinArray, numCellPins);
	// printHyperwts(hyperwts, numhyper);
    // printVertexSize(vertexSize, numCellsAndPads);
    cout << "\nPlz read README file first !!" << endl;
    cout << "It takes much time to run a large circuit file" << endl;
    cout << "All info (Q matrix, Q size, dx, dy, coordinates, WL) will be printed in output/*.txt file" << endl;
    cout << "coordinates.txt is info of before spreading and coordinates2.txt is info of after spreading" << endl;
    cout << "\nYou can just test all *.txt results which are in GUI)AllOutputs' dir, by GUI.py and jpg will be stored in 'parser/output/GUIplot.jpg'" << endl;
    cout << "Please compile ('make all') first before using GUI" << endl;
	// Determine the size of the Q matrix
	qMatrixSize = determineQMatrixSize(hEdge_idxToFirstEntryInPinArray, numhyper + 1, numCells_noPads);
    //cout << "Size of Q matrix: " << qMatrixSize << endl;

	// Allocate memory for the QMatrix
    QMatrix = new double*[qMatrixSize];
    for (int i = 0; i < qMatrixSize; ++i) {
        QMatrix[i] = new double[qMatrixSize];
    }

    // Initialize the matrix with zeros
    for (int i = 0; i < qMatrixSize; ++i) {
        for (int j = 0; j < qMatrixSize; ++j) {
            QMatrix[i][j] = 0.0;
        }
    }
	//Fill the Q matrix
    fillQMatrix();
	//printQMatrix();

	//Print pin locations
	//printPinLocations();

	// Fill xd vector
    fillXdYd();
	//printXdYd();
	// Implemented Conjugate Gradient to get pre-spreading X and Y coordinates
	solutionX = conjugateGradient(xd);
	solutionY = conjugateGradient(yd);
	// cout << "\nAfter conjugate gradient: " << endl;
	//printXYCor();

    // Spreading part:
    // Initialize the bins
    initializeBins();
    // Perform cell shifting
    pair< vector<double>, vector<double> > newPositions = shiftCells();
    newXPositions = newPositions.first;
    newYPositions = newPositions.second;
	

	//*****run "python plot.py" or "python3 plot.py" if needs to see placement graph before spreading in output directory*****
	//*****run "python plot2.py" or "python3 plot2.py" if needs to see placement graph after spreading in output directory*****
    calculateTotalWidthLength();
	cout << "\nSqrt of sum of square WL (only b/w movable cells) before spreading: " << preWL << endl;
	cout << "\nSqrt of sum of square WL (only b/w movable cells) after spreading: " << afterWL << endl;
    //"output/coordinates.txt"
	printCoordinatesToFile();
    //"output/coordinates2.txt"
    printCoordinatesToFile2();

    free(pinLocations);
	free(hEdge_idxToFirstEntryInPinArray);
	free(cellPinArray);
	free(hyperwts);

	// Free the allocated memory of Q matrix
    for (int i = 0; i < qMatrixSize; ++i) {
        delete[] QMatrix[i];
    }
    delete[] QMatrix;
}
