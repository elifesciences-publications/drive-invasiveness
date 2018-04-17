## Overview ##

This repository contains data and all of the code required to reproduce the figures presented in: C Noble, B Adlam, GM Church, KM Esvelt, MA Nowak, Current CRISPR gene drives systems are likely to be highly invasive in wild populations, bioRxiv (2017).

The code included comes in two types:

1. Code for running numerical simulations and generating data (written in a mixture of C++ and Matlab)
2. Code for loading the numerical simulation data and plotting it to reproduce the figures as shown in the text (written entirely in Matlab)

The following description addresses these two steps separately.

## Generating numerical simulation data ##

For Figures 1B, 1D, 1E-F and 2C, the data is generated directly in Matlab by the functions described later for plotting, so no additional step is required here. For Figures 2D-E and Figures 3-9, the data is generated by C++ code, which can be found in the directory "Data_Generation/". The C++ code runs simulations and saves data specifically for each of these figures, and the data is stored in subdirectories named "Data_Storage/Figure_X", where each subdirectory is named according to the figure the data goes with (e.g., "Data_Storage/Figure_3/").

Data corresponding to all of the figures in the text can be downloaded directly from the "Data_Storage/" subdirectories here in the repository and plotted as descriebd in the next step without re-generating it via the C++ code. If you'd like to run the simulations yourself, however, directions are below.

To use this C++ code, follow the below steps:

1. Build the C++ project contained in the directory "Data_Generation/". This was originally done using the CLion IDE, version 2017.1.2 on Mac OSX El Capitan, version 10.11.6. Note that the source code might require some edits for other operating systems (e.g., making new directories in the C++ code).
2. All data required for each of the figures listed above is automatically generated by functions in main.cpp, with the relevant figure indicated by the function name, e.g., "generate_data_fig_4()". So, Running main() exactly as written will produce all of the data for all of the figures successively; note that this will likely take a long time (roughly 3 days or so). Alternatively, comment out lines in main() that correspond to figures that you're not looking to generate; the times roughly required for each figure are indicated by comments in main().
3. All data produced by the C++ code in this way is saved in automatically-generated subdirectories of the directory called "Data_Storage/". The subdirectories are named according to the relevant figure. For example, data for Figure 3 will be found in the directory "Data_Storage/Figure_3/".

## Plotting data from numerical simulations to reproduce the figures ##

Once the data is created for a figure that you're interested in generating, navigate to the directory "Figure_Generation/". There, Matlab code can be found for automatic generation of each figure, with the related figure indicated in the filename of the Matlab code. For example, the code for generating the left panel of Figure 7 is "gen_figure_7_left.m".

A few things to note:
1. The data for Figures 1B, 1D, 1E-F and 2C is generated directly via these plotting functions, so these plots will take longer to generate. (The other plots should be made more-or-less instantaneously, since the required data will have already been generated via the C++ code as  described above.
2. Figure 10 is generated by running the code to generate Figure 1D (gen_figure_1D.m), since Fig. 10 is simply Fig. 1D with deterministic simulations overlaid, which is included by default in the code that generates Fig. 1D.
3. Subdirectories "Figure_Generation/Plotting Utilities/" and "Figure_Generation/Simulation Functions/" must be on the Matlab path, since they contain code called by the figure-generating functions.

## Contact information ##

Charleston Noble\
charleston.noble@gmail.com
