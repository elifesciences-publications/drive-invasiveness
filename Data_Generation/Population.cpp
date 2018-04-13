//
//  Created by Charleston Noble on 5/24/17.
//

#include "Population.hpp"

Population::Population(std::mt19937 *engine) : m_engine(engine) {
    Parameters a;
    m_size = a.N;
    m_uniformDist = std::uniform_real_distribution<>(0, 1);
    m_maxFrequency = 0;
    GenInitialVec();
    UpdateStatus();
}

Population::Population(Parameters a, std::mt19937 *engine)
        : m_a(a), m_size(a.N), m_engine(engine) {
    m_uniformDist = std::uniform_real_distribution<>(0, 1);
    m_maxFrequency = 0;
    GenInitialVec();
    UpdateStatus();
}

Population::Population( const Population& other )
    : m_size(other.m_size),
      m_finished(other.m_finished),
      m_maxFrequency(other.m_maxFrequency),
      m_driveFrequency(other.m_driveFrequency),
      m_a(other.m_a),
      m_driveExtinct(other.m_driveExtinct),
      m_engine(other.m_engine)
{
    std::copy(other.m_vec,other.m_vec+6,m_vec);
}

void Population::GenInitialVec() {
    m_vec[0] = m_a.N - m_a.i - m_a.r;
    m_vec[1] = 0;
    m_vec[2] = 0;
    m_vec[3] = m_a.i;
    m_vec[4] = 0;
    m_vec[5] = m_a.r;
}

void Population::UpdateBS() {

    // This function performs a single population update if the number of
    // offspring per reproduction is specified to be greater than one.
    // Otherwise the function Population::Update() performs the updating.

    double x[6] = {0, 0, 0, 0, 0, 0};
    double tempPop[6];
    memcpy(tempPop,m_vec,6*sizeof(double));
    
    // Pick first mate
    EMult(m_vec, m_a.BirthRates, x);    // x is m_pop .* birth rates
    Normalize(x);                       // x is reproduction probs
    int ind1 = IdxFromProbs(x);
    tempPop[ind1]--;
    
    // Pick second mate
    EMult(tempPop, m_a.BirthRates, x);
    Normalize(x);
    int ind2 = IdxFromProbs(x);
    tempPop[ind2]--;
    
    // Calculate offspring
    int ofsp[m_a.bs];
    for (int i=0; i<m_a.bs; i++) {
        ofsp[i] = IdxFromProbs(m_a.MatingTbl[ind1][ind2]);
    }
    
    // Choose individual(s) for removal
    int rem[m_a.bs];
    for (int i=0; i<m_a.bs; i++) {
        memcpy(x,tempPop,6*sizeof(double));
        Normalize(x);
        rem[i] = IdxFromProbs(x);
        tempPop[rem[i]]--;
    }
    
    // Update the population vector
    for (int i=0; i<m_a.bs; i++) {
        m_vec[rem[i]]--;
        m_vec[ofsp[i]]++;
    }
    
    // Figure out whether we're done or should keep updating
    UpdateStatus();

}

void Population::UpdateInbreeding() {

    // This function performs a single population update if there is a
    // nonzero selfing probability. Otherwise the function Population::Update()
    // performs the updating.

    double x[6] = {0, 0, 0, 0, 0, 0};
    double tempPop[6];
    memcpy(tempPop,m_vec,6*sizeof(double));

    // Pick first mate
    EMult(m_vec, m_a.BirthRates, x);    // x is m_pop .* birth rates
    Normalize(x);                   // x is reproduction probs
    int ind1 = IdxFromProbs(x);
    tempPop[ind1]--;

    // With probability F, have first individual mate with itself
    double r = m_uniformDist(*m_engine);
    int ind2 = -1;
    if (r < m_a.F) {
        ind2 = ind1;
    } else {
        // Pick second mate
        EMult(tempPop, m_a.BirthRates, x);
        Normalize(x);
        ind2 = IdxFromProbs(x);
        tempPop[ind2]--;
    }

    // Calculate offspring
    int ofsp = IdxFromProbs(m_a.MatingTbl[ind1][ind2]);

    // Choose individual for removal
    Normalize(tempPop);
    int rem = IdxFromProbs(tempPop);
    //
    // Update the population vector
    m_vec[rem]--;
    m_vec[ofsp]++;

    // Figure out whether we're done or should keep updating
    UpdateStatus();

}

void Population::Update() {

    // This function performs a single reproduction event.
    
    double x[6] = {0, 0, 0, 0, 0, 0};
    double tempPop[6];
    memcpy(tempPop,m_vec,6*sizeof(double));
    
    // Pick first mate
    EMult(m_vec, m_a.BirthRates, x);    // x is m_pop .* birth rates
    Normalize(x);                       // x is reproduction probs
    int ind1 = IdxFromProbs(x);
    tempPop[ind1]--;
    
    // Pick second mate
    EMult(tempPop, m_a.BirthRates, x);
    Normalize(x);
    int ind2 = IdxFromProbs(x);
    tempPop[ind2]--;
    
    // Calculate offspring
    int ofsp = IdxFromProbs(m_a.MatingTbl[ind1][ind2]);
    
    // Choose individual for removal
    EMult(tempPop, m_a.DeathRates, x);    // x is tempPop .* death rates
    Normalize(x);                         // x is death probs
    int rem = IdxFromProbs(x);

    // Update the population vector
    m_vec[rem]--;
    m_vec[ofsp]++;
    
    // Figure out whether we're done or should keep updating
    UpdateStatus();
    
}

void Population::RemoveIndividual(int genotype){
    m_vec[genotype]--;
    UpdateStatus();
}

void Population::AddIndividual(int genotype) {
    m_vec[genotype]++;
    UpdateStatus();
}

double Population::GetSize() {
    UpdateSize();
    return m_size;
};

int Population::GetUniformlyChosenIndividual() {

    double arr[6];
    double cumsum = 0;
    for (int i=0; i<6; i++) {cumsum+=m_vec[i]; arr[i] = cumsum;};
    for (int i=0; i<6; i++) {arr[i]/=cumsum;};

    double r = m_uniformDist(*m_engine);
    int sum = 0;
    for (int i=0; i<6; i++) {
        if (r > arr[i]) {
            sum++;
        }
    }
    return sum;

}

double Population::GetMaxFrequency() {
    return m_maxFrequency;
}

double Population::GetDriveFrequency() {
    return m_driveFrequency;
}

void Population::UpdateDriveFrequency() {
    double dCount = 2*m_vec[3] + m_vec[1] + m_vec[4];
    m_driveFrequency = dCount / (2 * m_size);
}

void Population::UpdateMaxFrequency() {
    if (m_driveFrequency > m_maxFrequency) {
        m_maxFrequency = m_driveFrequency;
    }
}

void Population::EMult(double* x, double* y, double* z) {
    z[0] = x[0]*y[0];
    z[1] = x[1]*y[1];
    z[2] = x[2]*y[2];
    z[3] = x[3]*y[3];
    z[4] = x[4]*y[4];
    z[5] = x[5]*y[5];
}

void Population::Normalize(double* x) {
    double sum;
    sum = x[0]+x[1]+x[2]+x[3]+x[4]+x[5];
    x[0] = x[0] / sum;
    x[1] = x[1] / sum;
    x[2] = x[2] / sum;
    x[3] = x[3] / sum;
    x[4] = x[4] / sum;
    x[5] = x[5] / sum;
}

int Population::IdxFromProbs(double* x) {
    double cumSum[6];
    cumSum[0] = x[0];
    cumSum[1] = cumSum[0] + x[1];
    cumSum[2] = cumSum[1] + x[2];
    cumSum[3] = cumSum[2] + x[3];
    cumSum[4] = cumSum[3] + x[4];
    cumSum[5] = cumSum[4] + x[5];
    
    double r = m_uniformDist(*m_engine);
    int sum = 0;
    for (int i=0; i<6; i++) {
        if (r > cumSum[i]) {
            sum++;
        }
    }
    return sum;
}

bool Population::DriveIsExtinct() {
    return m_driveExtinct;
}

void Population::UpdateStatus() {
    UpdateSize();
    int N = (int)m_size;
    m_finished = ( m_vec[0] == N || m_vec[1] == N || m_vec[2] == N ||
                   m_vec[3] == N || m_vec[4] == N || m_vec[5] == N );
    m_driveExtinct = ( m_vec[1] == 0 && m_vec[3] == 0 && m_vec[4] == 0 );
    UpdateDriveFrequency();
    UpdateMaxFrequency();
}

bool Population::IsFinished() {
    return m_finished;
}

double Population::CalcTotalFitness() {

    double total_fitness = 0;
    total_fitness += m_vec[0] * m_a.BirthRates[0];
    total_fitness += m_vec[1] * m_a.BirthRates[1];
    total_fitness += m_vec[2] * m_a.BirthRates[2];
    total_fitness += m_vec[3] * m_a.BirthRates[3];
    total_fitness += m_vec[4] * m_a.BirthRates[4];
    total_fitness += m_vec[5] * m_a.BirthRates[5];
    return total_fitness;

}

void Population::UpdateSize() {
    double total = 0;
    total += m_vec[0];
    total += m_vec[1];
    total += m_vec[2];
    total += m_vec[3];
    total += m_vec[4];
    total += m_vec[5];
    m_size = total;
}

bool Population::AllIndividualsSingleGenotype() {
    int types = 0;
    for (int i=0; i<6; i++) {
        if (m_vec[i] > 0) {
            types++;
        }
    }
    return (types == 1);
}

double Population::GetHeterozygoteCount() {
    return m_vec[1] + m_vec[2] + m_vec[4];
}