//
//  Created by Charleston Noble on 5/24/17.
//

#ifndef Population_hpp
#define Population_hpp

#include <random>
#include "Parameters.hpp"

class Population {
    
public:

    // Constructors
    Population(std::mt19937 *engine);
    Population(Parameters a, std::mt19937 *engine);
    Population(const Population & other);

    // Public member objects
    Parameters      m_a;

    // Member functions
    void        GenInitialVec();
    void        UpdateBS();
    void        Update();
    void        UpdateInbreeding();
    double      CalcTotalFitness();
    void        UpdateSize();
    void        RemoveIndividual(int genotype);
    void        AddIndividual(int genotype);
    double      GetSize();
    double      GetMaxFrequency();
    double      GetDriveFrequency();
    void        UpdateDriveFrequency();
    bool        DriveIsExtinct();
    int         GetUniformlyChosenIndividual();
    bool        AllIndividualsSingleGenotype();
    double      GetHeterozygoteCount();
    bool        IsFinished();

private:

    // Member objects
    std::mt19937    *m_engine;
    std::uniform_real_distribution<>    m_uniformDist;
    double          m_size;
    double          m_vec[6];
    bool            m_finished;
    double          m_maxFrequency;
    double          m_driveFrequency;
    bool            m_driveExtinct;

    // Member functions
    void    EMult(double x[6], double y[6], double* z);
    void    Normalize(double x[6]);
    int     IdxFromProbs(double* x);
    void    UpdateStatus();
    void    UpdateMaxFrequency();
    
};

#endif /* Population_hpp */