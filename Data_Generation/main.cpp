//
//  Created by Charleston Noble on 6/15/17.
//

#include "Simulation.hpp"
#include "StructPopSim.h"

void generate_data_fig_2D(std::string folder);
void generate_data_fig_2E(std::string folder);
void generate_data_fig_3(std::string folder);
void generate_data_fig_4(std::string folder);
void generate_data_fig_5(std::string folder);
void generate_data_fig_6(std::string folder);
void generate_data_fig_7_left(std::string folder);
void generate_data_fig_7_right(std::string folder);
void generate_data_fig_8(std::string folder);
void generate_data_fig_9(std::string folder);
void CreateDirectory(std::string folder);

int main(int argc, char *argv[]) {

    // All data will be saved in subdirectories of the directory written below.
    // The Matlab code for plotting the data assumes that it's found here, so
    // it's recommended to keep this unchanged.
    std::string directory = "../../Data_Storage/";

    // The directory will be created if it doesn't exist. Note, this code was
    // written for Mac, so the CreateDirectory function might not work on other OS's.
    CreateDirectory(directory);

    // Each function below runs code which generates data for the indicated figures,
    // numbered in the order that they're found in the text. Comment lines out if
    // you'd like to generate data for only some particular figures.
    generate_data_fig_2D(directory);
    generate_data_fig_2E(directory);
    generate_data_fig_3(directory);
    generate_data_fig_4(directory);
    generate_data_fig_5(directory);
    generate_data_fig_6(directory);
    generate_data_fig_7_left(directory);
    generate_data_fig_7_right(directory);
    generate_data_fig_8(directory);
    generate_data_fig_9(directory);

    // Return
    return 0;

}

void generate_data_fig_2D(std::string folder) {

    std::cout << "Starting calculations for Figure 2D..." << std::endl;

    folder += "Figure_2D/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    for (int i=0; i<51; i++) {
        StructPopSim::EscapeProbability(i, &engine, folder);
        std::cout << "Done with " << i + 1 << " out of " << 51 <<
                  " homing efficiency values." << std::endl;
    }
}

void generate_data_fig_2E(std::string folder) {

    std::cout << "Starting calculations for Figure 2E..." << std::endl;

    folder += "Figure_2E/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    for (int i=0; i<51; i++) {
        StructPopSim::InvasionProbability(i, &engine, folder);
        std::cout << "Done with " << i + 1 << " out of " << 51 <<
                  " migration rate values." << std::endl;
    }
}

void generate_data_fig_3(std::string folder) {

    std::cout << "Starting calculations for Figure 3..." << std::endl;

    folder += "Figure_3/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::NDependence(&engine, folder);
}

void generate_data_fig_4(std::string folder) {

    std::cout << "Starting calculations for Figure 4..." << std::endl;

    folder += "Figure_4/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::SGVDependence(&engine, folder);
}

void generate_data_fig_5(std::string folder) {

    std::cout << "Starting calculations for Figure 5..." << std::endl;

    folder += "Figure_5/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::FamilySizeDependenceReleaseSizeVaried(&engine, folder);
}

void generate_data_fig_6(std::string folder) {

    std::cout << "Starting calculations for Figure 6..." << std::endl;

    folder += "Figure_6/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::FamilySizeDependenceReleaseSizeConstant(&engine, folder);
}

void generate_data_fig_7_left(std::string folder) {

    std::cout << "Starting calculations for Figure 7 (left)..." << std::endl;

    folder += "Figure_7_left/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::fPDependence(&engine, folder);
}

void generate_data_fig_7_right(std::string folder) {

    std::cout << "Starting calculations for Figure 7 (right)..." << std::endl;

    folder += "Figure_7_right/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::fPDependenceDeathRate(&engine, folder);
}

void generate_data_fig_8(std::string folder) {

    std::cout << "Starting calculations for Figure 8..." << std::endl;

    folder += "Figure_8/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::fsDependence(&engine, folder);
}

void generate_data_fig_9(std::string folder) {

    std::cout << "Starting calculations for Figure 9..." << std::endl;

    folder += "Figure_9/";
    CreateDirectory(folder);

    std::random_device rd;
    std::mt19937 engine(rd());
    Simulation::IBDependence(&engine, folder);
}

void CreateDirectory(std::string folder)
{
    const int dir = system(std::string(std::string("mkdir ") + folder).c_str());
    if (dir < 0) {
        printf("Error creating directory.");
        exit(1);
    }
}