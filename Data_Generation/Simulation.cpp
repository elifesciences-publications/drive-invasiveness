//
//  Created by Charleston Noble on 5/24/17.
//

#include <iostream>
#include <array>
#include <sstream>
#include <cmath>
#include <fstream>
#include "Population.hpp"
#include "Simulation.hpp"
#include "StructPopSim.h"


Simulation::Simulation(Parameters a_in, std::mt19937 *engine)
        : m_a(a_in),
          m_pop(Population(m_a, engine)),
          m_engine(engine){}

void Simulation::Run() {

    // This function runs a single simulation. Different update functions
    // are used if there is selfing or if there are multiple offspring
    // per reproduction event, hence the conditions below.

    if (m_pop.m_a.bs == 1) {
        if (m_pop.m_a.F > 0) {
            while (!m_pop.IsFinished()) {
                m_pop.UpdateInbreeding();
            }
        } else {
            while (!m_pop.IsFinished()) {
                m_pop.Update();
            }
        }
    } else {
        if (m_pop.m_a.F == 0) {
            while (!m_pop.IsFinished()) {
                m_pop.UpdateBS();
            }
        } else {
            throw "Can't have inbreeding and offspring # greater than 1.";
        }
    }
}

void Simulation::NDependence(std::mt19937 *engine, std::string folder) {

    // This function generates the data for Figure 3.

    // Specify parameters for the figure.
    int simulations = 1000;
    double q = 1;
    double P = 0.5;
    double c = 0.1;
    double s = 0;
    int N[5] = {500, 1000, 2500, 5000, 10000};
    double i_arr[10] = {1,2,3,4,5,6,7,8,9,10};
    double r = 0;
    int bs = 1;
    double F = 0;
    bool birthRateUpdating = 0;

    // Array to save the maximum frequency distribution
    double maxFreqDist[simulations];

    int counter = 0;

    // Loop over each of the population sizes
    for(int N_i=0; N_i<sizeof(N)/sizeof(N[0]); N_i++) {

        // Get the population size we're currently testing
        int N_val = N[N_i];

        // Loop over each of the initial release sizes
        for(int i_idx=0; i_idx<sizeof(i_arr)/sizeof(i_arr[0]); i_idx++) {

            // Make a Parameters object based on the parameters above
            Parameters a(q, P, c, s, N_val, i_arr[i_idx], r, bs, F, birthRateUpdating);

            // Run a bunch of simulations and save the maximum frequency from each
            for(int s_i=0; s_i<simulations; s_i++) {
                Simulation sim(a, engine);                      // Initialize simulation
                sim.Run();                                      // Run simulation
                maxFreqDist[s_i] = sim.m_pop.GetMaxFrequency(); // Save maximum frequency
            }

            // Save to file
            std::ostringstream oss;
            oss << folder << "data_N_" << N[N_i] << "_i_" << i_arr[i_idx] << ".csv";
            SaveDistribution(oss.str(), maxFreqDist, simulations);

            // Display progress
            counter++;
            std::cout << "Done with " << counter << " out of " <<
                      sizeof(N)/sizeof(N[0]) * sizeof(i_arr)/sizeof(i_arr[0]) << " distributions.";
            std::cout << std::endl;
        }
    }
}

void Simulation::SGVDependence(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 4.

    // Specify parameters for the figure
    int simulations = 5000;
    double q = 1;
    double P = 0.5;
    double c = 0.1;
    double s = 0;
    double N = 1000;
    double i = 15;
    int bs = 1;
    double r_arr[11] = {0,50,100,150,200,250,300,350,400,450,500};
    double F = 0;

    // Array to save the maximum frequency distribution
    double maxFreqDist[simulations];

    // Loop over each of the SGV values
    for(int r_i=0; r_i<sizeof(r_arr)/sizeof(r_arr[0]); r_i++) {

        // Create a Parameters object based on the parameters above
        Parameters a(q, P, c, s, N, i, r_arr[r_i], bs, F);

        // Run a bunch of simulations
        for(int s_i=0; s_i<simulations; s_i++) {
            Simulation sim(a, engine);                      // Initialize simulation
            sim.Run();                                      // Run simulation
            maxFreqDist[s_i] = sim.m_pop.GetMaxFrequency(); // Save max frequency
        }

        // Save the distribution to a file
        std::ostringstream oss;
        oss << folder << "data_r_" << (double)r_arr[r_i]/N << ".csv";
        SaveDistribution(oss.str(), maxFreqDist, simulations);

        // Display progress
        std::cout << "Done with " << r_i + 1 << " out of " << 11
                  << " initial resistance frequencies." << std::endl;
    }
}

void Simulation::FamilySizeDependenceReleaseSizeVaried(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 5.

    // Specify parameters for the figure
    int simulations = 5000;
    double q = 1;
    double P = 0.5;
    double c = 0.1;
    double s = 0;
    double N[25]; for(int k=0; k<25; k++) {N[k]=250.0*(2.0*(double)(k+1)+6.0)/4.0;};
    double i[25]; for(int k=0; k<25; k++) {i[k]=8.0  *(2.0*(double)(k+1)+6.0)/4.0;};
    int bs_arr[25]; for(int k=0; k<25; k++) {bs_arr[k]=(k+1);};
    double r = 0;
    double F = 0;

    // Array to save the maximum frequency distribution
    double maxFreqDist[simulations];

    // Iterate over all the offspring counts (i.e., values of "k" from the paper)
    for(int bs_i=0; bs_i<sizeof(bs_arr)/sizeof(bs_arr[0]); bs_i++) {

        // Generate Parameters object based on the parameters above
        Parameters a(q, P, c, s, N[bs_i], i[bs_i], r, bs_arr[bs_i], F);

        // Run a bunch of simulations
        for(int s_i=0; s_i<simulations; s_i++) {
            Simulation sim(a, engine);                      // Initialize
            sim.Run();                                      // Run
            maxFreqDist[s_i] = sim.m_pop.GetMaxFrequency(); // Save max frequency
        }

        // Save the distribution to a file
        std::ostringstream oss;
        oss << folder << "data_k_" << (double)bs_arr[bs_i] << ".csv";
        SaveDistribution(oss.str(), maxFreqDist, simulations);

        // Display progress
        std::cout << "Done with " << bs_i+1 << " out of " << sizeof(bs_arr)/sizeof(bs_arr[0])
                  << " values of k." << std::endl;
    }
}

void Simulation::FamilySizeDependenceReleaseSizeConstant(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 6.

    // Specify parameters for the figure
    int simulations = 5000;
    double q = 1;
    double P = 0.5;
    double c = 0.1;
    double s = 0;
    double N = 500;
    double i = 15;
    int bs_arr[25]; for(int k=0; k<25; k++) {bs_arr[k]=(k+1);};
    double r = 0;
    double F = 0;

    // Array to save the maximum frequency distribution
    double maxFreqDist[simulations];

    // Loop over all the numbers of offspring per reproductive event ("k" from the paper)
    for(int bs_i=0; bs_i<sizeof(bs_arr)/sizeof(bs_arr[0]); bs_i++) {

        // Generate a Parameters object based on the parameters from above
        Parameters a(q, P, c, s, N, i, r, bs_arr[bs_i], F);

        // Run a bunch of simulations
        for(int s_i=0; s_i<simulations; s_i++) {
            Simulation sim(a, engine);                      // Initialize
            sim.Run();                                      // Run
            maxFreqDist[s_i] = sim.m_pop.GetMaxFrequency(); // Save max frequency
        }

        // Save the distribution to a file
        std::ostringstream oss;
        oss << folder << "data_k_" << (double)bs_arr[bs_i] << ".csv";
        SaveDistribution(oss.str(), maxFreqDist, simulations);

        // Display progress
        std::cout << "Done with " << bs_i+1 << " out of " << sizeof(bs_arr)/sizeof(bs_arr[0])
                  << " values of k." << std::endl;
    }
}

void Simulation::fPDependence(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 7 (left).

    // Specify parameters for the figure
    int simulations = 100;
    std::string fname = folder + "data.csv";
    double q = 1;
    const int dim1 = 51;
    const int dim2 = 51;
    double P_arr[dim1]; for(int i=0; i<dim1; i++) {P_arr[i] = (double)(i * 1.0/(dim1-1));}
    double c_arr[dim2]; for(int i=0; i<dim2; i++) {c_arr[i] = (double)(i * 0.5/(dim2-1));}
    double s = 0;
    double N = 500;
    double i0 = 15;
    int bs = 1;
    double r = 0;
    double F = 0;

    // Array to save the maximum frequency values
    double **maxFreqArr = new double *[dim1]; for(int i=0; i<dim1; i++) {maxFreqArr[i]=new double[dim2];};
    for (int i=0; i<dim1; i++) {for (int j=0; j<dim2; j++) {maxFreqArr[i][j]=-1;}};

    // Loop over the homing efficiency values (P)
    for (int i=0; i<dim1; i++) {
        double P = P_arr[i];

        // Loop over the drive fitness cost values (c)
        for (int j=0; j<dim2; j++) {
            double c = c_arr[j];

            // Generate a Parameters object based on the parameters from above
            Parameters a(q, P, c, s, N, i0, r, bs, F);

            // Get the average max frequency by totaling and dividing by simulations
            double total = 0;
            for (int k=0; k<simulations; k++) {
                Simulation sim(a, engine);              // Initialize simulation
                sim.Run();                              // Run simulation
                total += sim.m_pop.GetMaxFrequency();   // Count max frequency
            }
            maxFreqArr[i][j] = total / simulations;

            // Save the max frequency array after every value
            SaveArray(fname, maxFreqArr, dim1, dim2);
        }

        // Display progress
        std::cout << "Done with " << i+1 << " out of " << dim1 << " homing efficiency values."
                  << std::endl;
    }
}

void Simulation::fPDependenceDeathRate(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 7 (right).

    // Specify parameters for the figure
    int simulations = 100;
    std::string fname = folder + "data.csv";
    double q = 1;
    const int dim1 = 51;
    const int dim2 = 51;
    double P_arr[dim1]; for(int i=0; i<dim1; i++) {P_arr[i] = (double)(i * 1.0/(dim1-1));}
    double c_arr[dim2]; for(int i=0; i<dim2; i++) {c_arr[i] = (double)(i * 0.5/(dim2-1));}
    double s = 0;
    double N = 500;
    double i0 = 15;
    int bs = 1;
    double r = 0;
    double F = 0;
    bool birthRateUpdating = 0;

    // Array to save the maximum frequency values
    double **maxFreqArr = new double *[dim1]; for(int i=0; i<dim1; i++) {maxFreqArr[i]=new double[dim2];};
    for (int i=0; i<dim1; i++) {for (int j=0; j<dim2; j++) {maxFreqArr[i][j]=-1;}};

    // Loop over the homing efficiency (P) values
    for (int i=0; i<dim1; i++) {
        double P = P_arr[i];

        // Loop over the drive fitness cost (c) values
        for (int j=0; j<dim2; j++) {
            double c = c_arr[j];

            // Generate a Parameters object based on the parameters from above
            Parameters a(q, P, c, s, N, i0, r, bs, F, birthRateUpdating);

            // Get the average max frequency by totaling and dividing by simulations
            double total = 0;
            for (int k=0; k<simulations; k++) {
                Simulation sim(a, engine);              // Initialize simulation
                sim.Run();                              // Run
                total += sim.m_pop.GetMaxFrequency();   // Count max frequency
            }
            maxFreqArr[i][j] = total / simulations;

            // Save after every value is calculated
            SaveArray(fname, maxFreqArr, dim1, dim2);
        }

        // Display progress
        std::cout << "Done with " << i+1 << " out of " << dim1 << " homing efficiency values." << std::endl;
    }
}

void Simulation::fsDependence(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 8.

    // Specify parameters for the figure
    int simulations = 100;
    std::string fname = folder + "data.csv";
    double q = 1;
    const int dim1 = 51;
    const int dim2 = 51;
    double P = 0.5;
    double s_arr[dim1]; for(int i=0; i<dim1; i++) {s_arr[i] = (double)(i * 1.0/(dim1-1));}
    double c_arr[dim2]; for(int i=0; i<dim2; i++) {c_arr[i] = (double)(i * 0.5/(dim2-1));}
    double N = 500;
    double i0 = 15;
    int bs = 1;
    double r = 0;
    double F = 0;

    // Array to save the maximum frequency values
    double **maxFreqArr = new double *[dim1]; for(int i=0; i<dim1; i++) {maxFreqArr[i]=new double[dim2];};
    for (int i=0; i<dim1; i++) {for (int j=0; j<dim2; j++) {maxFreqArr[i][j]=-1;}};

    // Loop over the resistance fitness cost (s) values
    for (int i=0; i<dim1; i++) {
        double s = s_arr[i];

        // Loop over the drive fitness cost (c) values
        for (int j=0; j<dim2; j++) {
            double c = c_arr[j];

            // Generate a Parameters object using the parameters from above
            Parameters a(q, P, c, s, N, i0, r, bs, F);

            // Calculate the average max frequency by summing over simulations and dividing by simulations
            double total = 0;
            for (int k=0; k<simulations; k++) {
                Simulation sim(a, engine);              // Initialize simulation
                sim.Run();                              // Run
                total += sim.m_pop.GetMaxFrequency();   // Count max frequency
            }
            maxFreqArr[i][j] = total / simulations;

            // Save to file after every value is calculated
            SaveArray(fname, maxFreqArr, dim1, dim2);
        }

        // Display progress
        std::cout << "Done with " << i+1 << " out of " << dim1 << " resistance cost values." << std::endl;
    }
}

void Simulation::IBDependence(std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 9.

    // Specify parameters for the figure
    int simulations = 1000;
    double q = 1;
    const int dim1 = 3;
    const int dim2 = 51;
    double P_arr[3] = {0.15, 0.5, 0.9};
    double c = 0.1;
    double s = 0;
    double N = 500;
    double i0 = 15;
    int bs = 1;
    double r = 0;
    double F_arr[dim2]; for(double i=0; i<dim2; i++) {F_arr[(int)i] = (double)(i * 1/(dim2-1));}

    // Array to save the maximum frequency distribution
    double maxFreqDist[simulations];

    // Loop over the homing efficiency (P) values
    for (int i=0; i<dim1; i++) {
        double P = P_arr[i];

        // Loop over the selfing probability (F, here) values
        for (int j=0; j<dim2; j++) {
            double F = F_arr[j];

            // Generate a Parameters object based on the parameters above
            Parameters a(q, P, c, s, N, i0, r, bs, F);

            // Run a bunch of simulations
            for (int k=0; k<simulations; k++) {
                Simulation sim(a, engine);                      // Initialize
                sim.Run();                                      // Run
                maxFreqDist[k] = sim.m_pop.GetMaxFrequency();   // Save maximum frequency
            }

            // Save the distribution to a file
            std::ostringstream oss;
            oss << folder << "data_P_" << P << "_F_" << F << ".csv";
            SaveDistribution(oss.str(), maxFreqDist, simulations);
        }

        // Display progress
        std::cout << "Done with " << i+1 << " out of " << dim1 << " homing efficiency values." << std::endl;
    }
}

void Simulation::SaveArray(std::string fname, double **arr, int dim1, int dim2) {

    // This function saves a 2D array to a csv file.

    std::ofstream f;
    f.open(fname);
    for (int i=0; i<dim1; i++) {
        for (int j=0; j<dim2; j++) {
            f << arr[i][j] << ",";
        }
        f << "\n";
    }
    f.close();
}

void Simulation::SaveDistribution(std::string fname, double dist[], int sims) {

    // This function saves a 1D distribution to a csv file.

    std::ofstream f;
    f.open(fname);
    for (int i=0; i<sims; i++) {
        f << dist[i] << "\n";
    }
    f.close();
}