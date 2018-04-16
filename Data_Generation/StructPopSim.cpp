//
// Created by Charleston Noble on 6/15/17.
//

#include "StructPopSim.h"

StructPopSim::StructPopSim(Parameters a, int subpops, int subpop_size, double mig_rate,
                           std::mt19937 *engine, int invasion_cutoff_count)
        : m_pop(StructPop(subpops, subpop_size, (int)a.i, mig_rate, a, engine)),
          m_a(a),
          m_invasion_cutoff_count(invasion_cutoff_count),
          m_engine(engine) {}

StructPopSim::StructPopSim(Parameters a_in, int subpops_in, int subpop_size_in, double mig_rate,
                           std::mt19937 *engine)
        : m_pop(StructPop(subpops_in, subpop_size_in, (int)a_in.i, mig_rate, a_in, engine)),
          m_a(a_in),
          m_engine(engine),
          m_invasion_cutoff_count(1) {}

void StructPopSim::Run() {

    // This function runs a single simulation of the structured population model.

    // GetEscapeCount() includes only populations other than the initial release population
    while (m_pop.GetEscapeCount() < m_invasion_cutoff_count && !m_pop.m_drive_extinct) {
        m_pop.Update();
    }

    if (m_pop.m_invaded_pop_count >= m_invasion_cutoff_count) {
        m_escape_occurred = true;
    }

}

void StructPopSim::InvasionProbability(int index, std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 2E.

    // Specify the parameters for the figure
    int simulations = 10000;
    double q = 1;
    const int dim1 = 4;
    const int dim2 = 51;
    double c = 0.1;
    double s = 0;
    int i0 = 15;
    int bs = 1;
    double r = 0;
    double F = 0;
    double P = 0.5;
    int subpops = 5;
    int subpop_size = 100;

    // Array holding the number of populations invaded before "success"
    int invasion_count_arr[dim1]; for (int i=0; i<dim1; i++) {invasion_count_arr[i] = i+1;};

    // Filename for saving later
    std::ostringstream oss;
    oss << folder << "data_index_" << index << ".csv";

    // Migration rate vector
    double start = -5;
    double stop = log10(0.5);
    std::vector<double> mig_vals;
    generate_n(std::back_inserter(mig_vals), dim2, Logspace<>(start,stop,dim2));
    double m = mig_vals.at((size_t)index);

    // Array for saving the data
    double **escapeProbArr = new double *[dim1]; for(int i=0; i<dim1; i++) {escapeProbArr[i]=new double[1];};
    for (int j=0; j<dim1; j++) {escapeProbArr[j][0] = 0;};

    // Loop over the number of populations we're considering (i.e., the four different curves
    // seen in the figure.
    for (size_t j=0; j<dim1; j++) {

        // Get the number of populations required for "success"
        int invasion_count = invasion_count_arr[j];

        // Keep track of how many escapes happen over the simulations
        double escape_count = 0;

        // Run a bunch of simulations
        for (int sim=0; sim<simulations; sim++) {

            // Make a Parameters object based on the parameters above
            Parameters a(q, P, c, s, subpops*subpop_size, i0, r, bs, F);

            // Initialize and run a structured population simulation
            StructPopSim spsim(a, subpops, subpop_size, m, engine, invasion_count);
            spsim.Run();

            // Increment if an escape occurred
            if (spsim.m_escape_occurred) {
                escape_count++;
            }
        }

        // Define the escape probability as the number observed / simulations
        escapeProbArr[j][0] = escape_count / simulations;

        // Save the results
        SaveArray(oss.str(), escapeProbArr, dim1, 1);
    }
}


void StructPopSim::EscapeProbability(int index, std::mt19937 *engine, std::string folder) {

    // This function generates data for Figure 2D.

    // Specify the parameters for the figure
    int simulations = 1000;
    double q = 1;
    const int dim1 = 51;
    const int dim2 = 51;
    double c = 0.1;
    double s = 0;
    int i0 = 15;
    int bs = 1;
    double r = 0;
    double F = 0;
    int subpops = 5;
    int subpop_size = 100;

    // Homing efficiency array
    double P_arr[dim1]; for (int i=0; i<dim1; i++) {P_arr[i] = (double)i/(dim1-1);};
    double P = P_arr[index];

    // Filename for saving the data
    std::ostringstream oss;
    oss << folder << "data_index_" << index << ".csv";

    // Migration rate vector
    double start = -5;
    double stop = log10(0.5);
    std::vector<double> mig_vals;
    generate_n(std::back_inserter(mig_vals), dim2, Logspace<>(start,stop,dim2));

    // Array to save the data
    double **escapeProbArr = new double *[1]; for(int i=0; i<1; i++) {escapeProbArr[i]=new double[dim2];};
    for (int j=0; j<dim2; j++) {escapeProbArr[0][j] = 0;};

    for (size_t j=0; j<dim2; j++) {

        // Get the migration rate
        double m = mig_vals.at(j);

        // Count how many escapes occur over the various simulations
        double escape_count = 0;

        // Run a bunch of simulations
        for (int sim=0; sim<simulations; sim++) {

            // Generate a Parameters object based on the parameters above
            Parameters a(q, P, c, s, subpops*subpop_size, i0, r, bs, F);

            // Initialize and run a structured population simulation
            StructPopSim spsim(a, subpops, subpop_size, m, engine);
            spsim.Run();

            // If an escape occurred, count it
            if (spsim.m_escape_occurred) {
                escape_count++;
            }
        }

        // Escape probability is the number of escapes / simulations
        escapeProbArr[0][j] = escape_count / (double)simulations;

        // Save the results
        SaveArray(oss.str(), escapeProbArr, 1, dim2);

        // Display progress
        std::cout << "Done with " << j+1 << " out of " << dim2 << " migration rate values." << std::endl;
    }
}

void StructPopSim::SaveArray(std::string fname, double **arr, int dim1, int dim2) {

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