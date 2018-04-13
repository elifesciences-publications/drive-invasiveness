//
//  Created by Charleston Noble on 5/24/17.
//

#include "Parameters.hpp"
#include <iostream>

Parameters::Parameters() {
    q = 1;
    P = 0.5;
    c = 0.1;
    s = 0;
    N = 500;
    i = 10;
    r = 0;
    bs= 1;
    F = 0;
    birthRateUpdating = 1;

    GenGameteTable();
    GenMatingTable();
    SetBirthRates();
    SetDeathRates();
}

Parameters::Parameters(double q_in, double P_in, double c_in, double s_in,
                       double N_in, double i_in, double r_in, int bs_in, double F_in) {
    q = q_in;
    P = P_in;
    c = c_in;
    s = s_in;
    N = N_in;
    i = i_in;
    r = r_in;
    bs = bs_in;
    F = F_in;
    birthRateUpdating = 1;
    GenGameteTable();
    GenMatingTable();
    SetBirthRates();
    SetDeathRates();
}

Parameters::Parameters(double q_in, double P_in, double c_in, double s_in,
                       double N_in, double i_in, double r_in, int bs_in, double F_in,
                       bool birthRateUpdating_in) {
    q = q_in;
    P = P_in;
    c = c_in;
    s = s_in;
    N = N_in;
    i = i_in;
    r = r_in;
    bs = bs_in;
    F = F_in;
    birthRateUpdating = birthRateUpdating_in;
    GenGameteTable();
    GenMatingTable();
    SetBirthRates();
    SetDeathRates();
}

void Parameters::GenMatingTable() {

    // Produces a 3D array (6x6x6) whose (i,j,k) value is the probability that
    // parents having types i and j produce an offspring of type k. The order of
    // types is WW, WD, WR, DD, DR, RR.
    
    int ta[6][2] = {{0, 0},{0, 1},{0, 2},{1, 1},{1, 2},{2,2}};
    int t3a0;
    int t3a1;
    double prob;
    double prob_lr;
    double prob_rl;
    
    for (int i=0; i<=5; i++) {
        int type1 = i;
        
        for (int j=0; j<=5; j++) {
            int type2 = j;
            
            for (int k=0; k<=5; k++) {
                t3a0 = ta[k][0];
                t3a1 = ta[k][1];
                
                if (t3a0 == t3a1) {
                    prob = GameteTbl[type1][t3a0] * GameteTbl[type2][t3a1];
                } else {
                    prob_lr = GameteTbl[type1][t3a0] * GameteTbl[type2][t3a1];
                    prob_rl = GameteTbl[type1][t3a1] * GameteTbl[type2][t3a0];
                    prob = prob_lr + prob_rl;
                }
                MatingTbl[i][j][k] = prob;
            }
        }
    }
}

void Parameters::GenGameteTable() {

    // Produces a 6x3 array whose (i,j) entry is the probability that
    // an individual of type i produces a gamete of type j. Individuals
    // are ordered WW, WD, WR, DD, DR, RR, and gametes are ordered W, D, R.

    int t, g;
    for (t=0; t<=5; t++) {
        for (g=0; g<=2; g++) {
            if (g==0) {
                if (t==0) {
                    GameteTbl[t][g] = 1;
                }
                else if (t==1) {
                    GameteTbl[t][g] = (1-q)/2;
                }
                else if (t==2) {
                    GameteTbl[t][g] = 0.5;
                }
                else {
                    GameteTbl[t][g] = 0;
                }
            }
            else if (g == 1) {
                if (t==1) {
                    GameteTbl[t][g] = (1+q*P)/2;
                }
                else if (t==3) {
                    GameteTbl[t][g] = 1;
                }
                else if (t==4) {
                    GameteTbl[t][g] = 0.5;
                }
                else {
                    GameteTbl[t][g] = 0;
                }
            }
            else if (g == 2) {
                if (t==1) {
                    GameteTbl[t][g] = q*(1-P)/2;
                }
                else if (t==2) {
                    GameteTbl[t][g] = 0.5;
                }
                else if (t==4) {
                    GameteTbl[t][g] = 0.5;
                }
                else if (t==5) {
                    GameteTbl[t][g] = 1;
                }
                else {
                    GameteTbl[t][g] = 0;
                }
            }
        }
    }
}

void Parameters::SetBirthRates() {

    // Defines birth rates for the different genotypes, ordered
    // WW, WD, WR, DD, DR, RR. birthRateUpdating is a boolean value
    // which is 1 if matings happen proportional to fitness or 0
    // if deaths are varied instead.

    if (birthRateUpdating == 1) {
        BirthRates[0] = 1;
        BirthRates[1] = 1 - c;
        BirthRates[2] = 1;
        BirthRates[3] = 1 - c;
        BirthRates[4] = 1 - c;
        BirthRates[5] = 1 - s;
    } else {
        BirthRates[0] = 1;
        BirthRates[1] = 1;
        BirthRates[2] = 1;
        BirthRates[3] = 1;
        BirthRates[4] = 1;
        BirthRates[5] = 1;
    }

}

void Parameters::SetDeathRates() {

    // If birthRateUpdating is set to 1 then death rates are set
    // to be uniform, but if birthRateUpdating is 0 then death
    // rates are set to be the inverse of fitness.

    if (birthRateUpdating == 1) {
        DeathRates[0] = 1;
        DeathRates[1] = 1;
        DeathRates[2] = 1;
        DeathRates[3] = 1;
        DeathRates[4] = 1;
        DeathRates[5] = 1;
    } else {
        double use_c;
        if (c == 0) {
            use_c = 1.0e-6;
        } else {
            use_c = c;
        }
        double use_s;
        if (s == 0) {
            use_s = 1.0e-6;
        } else {
            use_c = s;
        }
        DeathRates[0] = 1;
        DeathRates[1] = 1/(1-c);
        DeathRates[2] = 1;
        DeathRates[3] = 1/(1-c);
        DeathRates[4] = 1/(1-c);
        DeathRates[5] = 1/(1-s);
    }
}