EE 5301 PA3  ANALYTIC PLACEMENT
Name: Yuan Jiang
ID: 5659433

Project Description:
    The project implements an analytic placement algorithm to partially solve the placement problem based on the FastPlace paper. It improves on traditional quadratic placement methods by reducing cell overlap and wirelength through three main strategies: Cell Shifting for efficient cell distribution, Iterative Local Refinement for optimizing wirelength, and a Hybrid Net Model to enhance solver speed by reducing non-zero entries in the connectivity matrix. The program will be given some large circuit files as a problem. After running it, all outputs will be store in txt files with width length and plot cells/IOPins placement in jpg file for before/after spreading. It also allows user to use GUI to plot by themselves to check the results.

Project results:
    Can't read ibm10 and ibm16 successfully.
    I checked my before spreading files of both of toy files, the Q matrix, Q size, dx, dy, x/y coordinates, WL and cell placement.jpg are all the same as Prof. Kia's result in forum. So I am pretty sure everything about Phase 1, Conjugate Gradient method and WL calculation algorithm are correct.
    For after spreading part in phase 2, I strictly followed paper's spreading algorithm and have checked carefully my funtion and implementations, the spreading results are reasonable.
    For ibm circuit files, the placement plots after spreading shows some cells are outside of I/O pads boundary, but I found it related to alpha value which is fixed to 0.8. When I toggled it to smaller value like 0.2, this issue doesn't happen. But in result plots, I keep placement plots with alpha = 0.8.

How to run:

    1. "make clean" -> "make all" -> "output" directory will be created
    2. "make test***" (*** is the circuit file name without extension) or "./suraj_parser ***"
    3. Then, "python3 plot.py" and "python3 plot2.py" to get plots of cell placements of before and after spreading correspondingly
    4. All info (Q matrix, Q size, dx, dy, coordinates, WL, cell placement.jpg) will be stored in "output" directory
    
    GUI:
    1. compile first (step 1 above)
    2. "python3 GUI.py" which will create a GUI window to select *.txt file to get plot of it in '/parser/output/GUIplot.jpg'

Comments:
    1. "AllWL.txt" contains all wirelength value for each benchmark before and after spreading
    2. "(GUI)AllOutputs" directory contains all benchmarks' coordinates files(.txt) and cell placement file(.jpg) before and after spreading.
    3. "(GUI)AllOutputs" directory is used for GUI.py in "parser" directory
    4. Before "python3 GUI.py", compile first and a "output" folder will be created
    5. After "python3 GUI.py", a GUI small window will be created and shown on screen, the user can select all "*.txt" in "(GUI)AllOutputs" directory to check plots, which will be plotted in '/parser/output/GUIplot.jpg'. The user can select multiple times of *.txt files and check GUIplot.jpg
    6. There are also .jpg cell placement results in "(GUI)AllOutputs"
    7. In "parser" directory, I changed Makefile to let "make all" create "output" directory and "make clean" delete "output"
    8. In "output" directory, after running the program, all info (Q matrix, Q size, dx, dy, coordinates, WL) will be stored in .txt file before and after spreading
    9. After get .txt files, run "python plot.py" and "python plot2.py" will create "plot.jpg" and "plot2.jpg" which are corresponding to cell placement of pre-spreading and after-spreading in "output"
    10. Certainly, the user can directly use GUI to check the results above without using plot*.py
    11. Terminal won't show anything except WL value of before spreading and after spreading