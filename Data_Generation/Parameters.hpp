//
//  Created by Charleston Noble on 5/24/17.
//

#ifndef Parameters_hpp
#define Parameters_hpp

class Parameters {
    
public:
    
    double MatingTbl[6][6][6];
    double GameteTbl[6][3];
    double BirthRates[6];
    double DeathRates[6];
    
    double P;   // Homing efficiency
    double q;   // Cutting efficiency
    double c;   // Drive fitness cost
    double s;   // Resistance fitness cost
    double N;   // m_pop size
    double i;   // Initial drive individuals
    double r;   // Initial resistant individuals
    int bs;     // "Brood size", number of offspring per reproduction (default 1)
    double F;   // Selfing probability
    bool birthRateUpdating;     // 1 if birth rate updating, 0 for death rate
    
    Parameters();
    Parameters(double q, double P, double c, double s, double N,
               double i, double r, int bs_in, double F_in);
    Parameters(double q, double P, double c, double s, double N,
               double i, double r, int bs_in, double F_in, bool bru_in);
    void GenMatingTable();
    void GenGameteTable();
    void SetBirthRates();
    void SetDeathRates();
    
};

#endif /* Parameters_hpp */