//
// Created by Charleston Noble on 6/15/17.
//

#ifndef MORANDRIVE_STRUCTPOPSIMULATION_H
#define MORANDRIVE_STRUCTPOPSIMULATION_H

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <iostream>
#include <array>
#include <string>
#include <sstream>
#include <fstream>
#include "Population.hpp"
#include "Simulation.hpp"
#include "StructPop.h"

class StructPopSim {

public:

    // Constructor
    StructPopSim(Parameters a, int subpops, int subpopsize, double mig_rate, std::mt19937 *engine);
    StructPopSim(Parameters a, int subpops, int subpop_size, double mig_rate,
                 std::mt19937 *engine, int invasion_cutoff_count);

    // Member objects
    int             m_invasion_cutoff_count;
    bool            m_escape_occurred = false;
    StructPop       m_pop;
    Parameters      m_a;
    std::mt19937    *m_engine;

    // Member functions
    void            Run();
    static void     EscapeProbability(int index, std::mt19937 *engine, std::string folder);
    static void     InvasionProbability(int index, std::mt19937 *engine, std::string folder);
    static void     SaveArray(std::string fname, double **arr, int dim1, int dim2);

};

// For generating log-spaced vectors
template<typename T = double>
class Logspace {
private:
    T curValue, base, step;
public:
    Logspace(T first, T last, int num, T base = 10.0) : curValue(first), base(base){
        step = (last - first)/(num-1);
    }
    T operator()() {
        T retval = pow(base, curValue);
        curValue += step;
        return retval;
    }
};

#endif //MORANDRIVE_STRUCTPOPSIMULATION_H