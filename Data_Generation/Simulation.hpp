//
//  Created by Charleston Noble on 5/25/17.
//

#ifndef Simulation_hpp
#define Simulation_hpp

#include <vector>
#include "Parameters.hpp"
#include "Population.hpp"
#include <array>


class Simulation {
    
public:

    // Constructor
    Simulation(Parameters a, std::mt19937 *engine);

    // Member objects
    std::vector<std::array<double,6>> m_allPops = {};

    // Member functions
    void        Run();
    static void NDependence(std::mt19937 *engine, std::string folder);
    static void SGVDependence(std::mt19937 *engine, std::string folder);
    static void FamilySizeDependenceReleaseSizeConstant(std::mt19937 *engine, std::string folder);
    static void FamilySizeDependenceReleaseSizeVaried(std::mt19937 *engine, std::string folder);
    static void fPDependence(std::mt19937 *engine, std::string folder);
    static void fPDependenceDeathRate(std::mt19937 *engine, std::string folder);
    static void fsDependence(std::mt19937 *engine, std::string folder);
    static void IBDependence(std::mt19937 *engine, std::string folder);
    static void SaveDistribution(std::string fname, double dist[], int sims);
    static void SaveArray(std::string fname, double **array, int dim1, int dim2);
    
private:

    // Member objects
    Parameters      m_a;
    Population      m_pop;
    std::mt19937    *m_engine;
    
};

#endif /* Simulation_hpp */