//
// Created by Charleston Noble on 6/15/17.
//

#ifndef MORANDRIVE_STRUCTUREDPOPULATION_H
#define MORANDRIVE_STRUCTUREDPOPULATION_H

#include "Population.hpp"
#include <vector>

class StructPop {

public:

    // Constructor
    StructPop(int subpops, int subpop_size, int initial_drives, double migration_rate,
              Parameters a, std::mt19937 *engine);

    // Member objects
    std::vector<Population>     m_whole_population;
    std::vector<double>         m_subpop_fitnesses;
    std::vector<double>         m_subpop_sizes;
    std::vector<int>            m_invaded_pops;
    Parameters                  m_params;
    int                         m_initial_drives;
    int                         m_subpops;
    int                         m_subpop_size;
    bool                        m_skip_reproductions;
    bool                        m_drive_extinct;
    double                      m_migration_rate;
    int                         m_invaded_pop_count;
    double                      m_invasion_frequency;
    int                         m_introduction_subpop;

    // Member functions
    void    Update();
    int     GetEscapeCount();

private:

    // Member objects
    std::mt19937 *m_engine;
    std::uniform_real_distribution<double> m_uniformDist;

    // Member functions
    void    MigrationEvent();
    void    ReproductionEvent();
    void    UpdateSubpopSizes();
    void    UpdateSubpopFitnesses();
    void    CheckEscape();
    void    CheckDriveExtinct();
    void    UpdateReproductionSkippingBool();
    int     GetSourcePopulation();
    int     GetDestPopulation(int source_pop);
    int     GetReproductionPop();
    void    UpdateStatus();

};


#endif //MORANDRIVE_STRUCTUREDPOPULATION_H
