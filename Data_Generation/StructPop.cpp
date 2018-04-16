//
// Created by Charleston Noble on 6/15/17.
//

#include "StructPop.h"
#include <iostream>

StructPop::StructPop(int subpops, int subpop_size, int initial_drives,
                     double m_mr, Parameters a, std::mt19937 *engine)
        : m_subpops(subpops),
          m_initial_drives(initial_drives),
          m_subpop_size(subpops),
          m_params(a),
          m_migration_rate(m_mr),
          m_invasion_frequency(0.1),
          m_engine(engine) {

    // Initialize the whole-population vector
    std::vector<Population> pop_vec;
    Parameters b = a;
    b.i = initial_drives;
    b.N = subpop_size;
    Population init_pop(b, engine);
    pop_vec.push_back(init_pop);
    m_introduction_subpop = 0;
    b.i = 0;
    for (int i=1; i<subpops; i++) {
        Population pop(b, engine);
        pop_vec.push_back(pop);
    }
    m_whole_population = pop_vec;

    // Initialize the subpopulation fitness vector
    std::vector<double> fitness_vec;
    for (int i=0; i<subpops; i++) {fitness_vec.push_back(0);};
    m_subpop_fitnesses = fitness_vec;
    UpdateSubpopFitnesses();

    // Initialize the subpopulation size vector
    std::vector<double> subpop_sizes;
    for (int i=0; i<subpops; i++) {subpop_sizes.push_back((double)subpop_size);};
    m_subpop_sizes = subpop_sizes;

    // Initialize the vector of invaded populations
    std::vector<int> invaded_pops;
    for (int i=0; i<subpops; i++) {
        invaded_pops.push_back(0);
    }
    m_invaded_pops = invaded_pops;
    m_invaded_pop_count = 0; // This will update in UpdateStatus();

    // Set the random number generator
    m_uniformDist = std::uniform_real_distribution<>(0, 1);

    // Set/update all the state variables
    UpdateStatus();
}

void StructPop::Update() {

    // This function performs a single step of the multiple-population process.

    double r = m_uniformDist(*m_engine);

    // Note that we skip reproduction events if the state of the population cannot
    // change via reproductions, since this drastically speeds up the calculations.
    // For example, an allele is fixed in each subpopultion, but no allele is
    // fixed globally. This is exactly correct for calculations where all we
    // care about are escape probabilities or maximum frequencies, etc. It's
    // incorrect if we care about the trajectory over time, so we instead do
    // those calculations in the Matlab version of the code.
    if (m_skip_reproductions) {
        MigrationEvent();
    }
    else {
        if (r < m_migration_rate) {
            MigrationEvent();
        } else {
            ReproductionEvent();
        }
    }

    // Update all the state variables
    UpdateStatus();

}

void StructPop::MigrationEvent() {

    // Pick the first population
    int subpop1 = GetSourcePopulation();

    // Pick the individual to move
    int genotype = m_whole_population.at((size_t)subpop1).GetUniformlyChosenIndividual();

    // Pick the second population
    int subpop2 = GetDestPopulation(subpop1);

    // Update the vectors
    m_whole_population.at((size_t)subpop1).RemoveIndividual(genotype);
    m_whole_population.at((size_t)subpop2).AddIndividual(genotype);

}

void StructPop::ReproductionEvent() {

    // Pick the population to do the reproduction in
    int subpop = GetReproductionPop();

    // Run an update step
    m_whole_population.at((size_t)subpop).Update();

}

int StructPop::GetEscapeCount() {
    CheckEscape();
    return m_invaded_pop_count;
}

void StructPop::CheckEscape() {

    // Check drive allele frequencies in all populations besides the initial-release population
    // Assumes the drive is initially released in population 0.

    for (size_t i=0; i<m_subpops; i++) {
        if (m_whole_population.at(i).GetMaxFrequency() > m_invasion_frequency) {
            if (m_invaded_pops.at(i) == 0) {
                m_invaded_pops.at(i) = 1;
                if (i != m_introduction_subpop) {
                    m_invaded_pop_count++;
                }
            }
        }
    }
}

void StructPop::CheckDriveExtinct() {

    // Checks whether the drive has gone completely extinct in all subpopulations.
    // Logic is that the drive is extinct if (it's extinct in population 1) and
    // (it's extinct in population 2) and (...), etc. So we start with true and
    // iterate over the subpopulations and check each.

    bool drive_extinct = true;
    for (size_t i=0; i<m_subpops; i++) {
        drive_extinct = drive_extinct && m_whole_population.at(i).DriveIsExtinct();
    }
    m_drive_extinct = drive_extinct;
}

int StructPop::GetDestPopulation(int source_pop) {

    // Samples a destination population for a migration (uniformly)
    // given a source population

    std::vector<double> cum_prob_vec;
    for (int i=0; i<m_subpops; i++) {
        cum_prob_vec.push_back((double)1/(m_subpops-1));
    };
    cum_prob_vec.at((size_t)source_pop) = 0;
    double cumsum = 0;
    for (size_t i=0; i<m_subpops; i++) {
            cumsum+=cum_prob_vec.at(i); cum_prob_vec.at(i)=cumsum;
    }

    double r = m_uniformDist(*m_engine);
    int sum = 0;
    for (size_t i=0; i<m_subpops; i++) {
        if (r > cum_prob_vec.at(i)) {
            sum++;
        }
    }
    return sum;

}

int StructPop::GetSourcePopulation() {

    // Samples a source population with probability proportional to its size

    std::vector<double> cum_size_vec;
    double cumsum = 0;
    for (size_t i=0; i<m_subpops; i++) {
        cumsum+=m_subpop_sizes.at(i);
        cum_size_vec.push_back(cumsum);
    };
    for (int i=0; i<m_subpops; i++) {
        cum_size_vec[i]/=cumsum;
    }

    double r = m_uniformDist(*m_engine);
    int sum = 0;
    for (size_t i=0; i<m_subpops; i++) {
        if (r > cum_size_vec.at(i)) {
            sum++;
        }
    }
    return sum;

}

int StructPop::GetReproductionPop() {

    // This function samples a population for reproduction with probability
    // proportion to its total fitness (more precisely, with probability
    // proportional to the square of the sum of its individuals' fitness).

    std::vector<double> cum_fit_vec;
    double cumsum = 0;
    for (size_t i=0; i<m_subpops; i++) {
        cumsum+=m_subpop_fitnesses.at(i);
        cum_fit_vec.push_back(cumsum);
    };
    for (int i=0; i<cum_fit_vec.size(); i++) {
        cum_fit_vec[i]/=cumsum;
    }

    double r = m_uniformDist(*m_engine);
    int sum = 0;
    for (size_t i=0; i<m_subpops; i++) {
        if (r > cum_fit_vec.at(i)) {
            sum++;
        }
    }
    return sum;
}

void StructPop::UpdateStatus() {

    // Updates all the state variables

    UpdateSubpopSizes();
    UpdateSubpopFitnesses();
    CheckEscape();
    CheckDriveExtinct();
    UpdateReproductionSkippingBool();

}

void StructPop::UpdateReproductionSkippingBool() {

    // Checks whether we should force migration events; this occurs iff
    // it is not possible for the state of the (whole) population to change
    // by mating alone--for example, if two different alleles are fixed
    // in different subpopulations. Then mating within subpopulations can't
    // change anything, yet the simulation isn't over because migration could
    // bring an allele to a population where it isn't fixed.

    // Check whether all individuals are single types in each subpop
    bool all_single_genotypes = true;
    for (size_t i=0; i<m_subpops; i++) {
        all_single_genotypes =
                all_single_genotypes &&
                m_whole_population.at(i).AllIndividualsSingleGenotype();
    }

    if (!all_single_genotypes) {
        m_skip_reproductions = false;
        return;
    }

    // Check whether all individuals are homozygotes
    double total_heterozygotes = 0;
    for (size_t i=0; i<m_subpops; i++) {
        total_heterozygotes += m_whole_population.at(i).GetHeterozygoteCount();
    }

    m_skip_reproductions = total_heterozygotes == 0;

}

void StructPop::UpdateSubpopSizes() {

    // Updates the sizes of the subpopulations. Used for calculating
    // probabilities of subpopulations serving as migration sources.

    for (size_t i=0; i<m_subpops; i++) {
        m_subpop_sizes.at(i) = m_whole_population.at(i).GetSize();
    }
}

void StructPop::UpdateSubpopFitnesses() {

    // Updates the total square of the sum of individuals' fitnesses
    // for each subpopulation

    for (size_t i=0; i<m_subpops; i++) {
        double temp = pow(m_whole_population.at(i).CalcTotalFitness(), 2);
        m_subpop_fitnesses.at(i) = temp;
    }
}