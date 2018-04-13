function [a, out] = moran_wm(i, P, sims, max_gens, output_type)

% Simulates the Moran process with sexual reproduction and three alleles: 
% WT, drive, and resistant. This function assumes one well-mixed
% population.
%
% In vectors used throughout, genotypes are ordered as follows:
%   1 - WW
%   2 - WD
%   3 - WR
%   4 - DD
%   5 - DR
%   6 - RR
%
% Gametes are ordered as follows:
%   1 - W
%   2 - D
%   3 - R
%
% output_type is a string which is either "full" or "max_dist". The former
% returns the entire population history, while the latter only keeps track
% of the maximum drive frequency observed in each simulation.

if strcmp(output_type, 'full')
    full_output = 1;
elseif strcmp(output_type, 'max_dist')
    full_output = 0;
else
    error('Invalid output type.')
end

% Generate a struct to store the parameters so they can be easily passed
% around.
a = struct();
a.q = 1.0;      % Cutting probability
a.c = 0.1;      % Fitness cost of drive
a.s = 0.0;      % Fitness cost of resistance
a.N = 500;      % Population size

if nargin == 0
    a.i = 15;       % Initial DD individuals
    a.P = 0.9;      % Homing efficiency
    a.sims = 1000;  % Number of simulations to run
    max_gens = 80;  % Maximum time to simulate for (in generations)
else
    a.i = i;
    a.P = P;
    a.sims = sims;
end
a.max_steps = a.N * max_gens;

tbl = mating_table(a);      % Mating table
birth = init_birth(a);      % Birth rates
death = init_death(a);      % Death rates

if full_output
    store_cell = cell(sims,1);
else
    store_dist = zeros(sims,1);
end

for sim = 1:sims
    
    % Initialize the population
    pop = init_pop(a);
    
    if full_output
        % Generate a vector to store the population state over time
        store_vec = zeros(a.max_steps,6);
        store_vec(1,:) = pop;
    else
        max_freq = 0;
    end
    len = 1;
    
    % Update the population until stopping condition reached
    while sum(pop == a.N) == 0 && len <= a.max_steps
        pop = update(pop,tbl,birth,death);
        if full_output
            store_vec(len+1,:) = pop;
        else
            d = drive_allele_count(pop);
            f = d / (2*a.N);
            if f > max_freq
                max_freq = f; 
            end
            if f == 0
                % If the drive has gone extinct, then stop tracking; max
                % frequency can't be higher than it currently is
                break
            end
        end
        len = len + 1;
    end
    
    if full_output
         % Store the allele counts over time in store_cell
        store_vec(len+1:end,:) = [];
        store_cell{sim} = allele_counts(store_vec);
    else
        % Store the maximum frequency
        store_dist(sim) = max_freq;
    end
end

if full_output
    out = store_cell;
else
    out = store_dist;
end

end


%=========================================================================%
function cts = allele_counts(pop)

% Converts genotype counts to allele counts.

cts = zeros(size(pop,1),2);
cts(:,1) = 2*pop(:,4) + pop(:,2) + pop(:,5);
cts(:,2) = 2*pop(:,6) + pop(:,3) + pop(:,5);

end


%=========================================================================%
function d = drive_allele_count(vec)

% Just returns the count of drive alleles given a population vector.

d = 2*vec(4) + vec(2) + vec(5);

end


%=========================================================================%
function pop_out = update(pop_in, tbl, birth, death)

% Runs one reproduction event and updates the population accordingly.

% Step 1: Pick two individuals to reproduce
rep_probs = (pop_in .* birth) / sum(pop_in .* birth);
ind1 = find(rand <= cumsum(rep_probs), 1);
temp = pop_in; temp(ind1) = temp(ind1) - 1;
rep_probs = (temp .* birth) / sum(temp .* birth);
ind2 = find(rand <= cumsum(rep_probs), 1);
temp(ind2) = temp(ind2) - 1;

% Step 2: Calculate offspring
offspring = find(rand <= cumsum(tbl(ind1,ind2,:)), 1);

% Step 3: Choose one individual for removal
rem_probs = (temp .* death) / sum(temp .* death);
removal = find(rand <= cumsum(rem_probs), 1);

% Step 4: Replace
pop_out = pop_in;
pop_out(removal) = pop_out(removal) - 1;
pop_out(offspring) = pop_out(offspring) + 1;

end


%=========================================================================%
function pop = init_pop(a)

% Initializes the population with i DD individuals and N-i WW individuals.

pop = zeros(6,1);
pop(4) = a.i;
pop(1) = a.N - a.i;

end


%=========================================================================%
function fit = init_death(~)

% Returns the death rates of each genotype in the order described at the 
% top of the file.

fit = ones(6,1);

end


%=========================================================================%
function fit = init_birth(a)

% Returns the birth rates of each genotype in the order described at the 
% top of the file.

fit = zeros(6,1);
fit(1) = 1;
fit(2) = 1-a.c;
fit(3) = 1;
fit(4) = 1-a.c;
fit(5) = 1-a.c;
fit(6) = 1-a.s;

end


%=========================================================================%
function tbl = mating_table(a)

% This function returns a 6x6x6 array, wherein the (i,j,k) entry is the
% probability of two parents of types i and j, respectively, producing an
% offspring of type k. This array should sum to 1 over the third dimension.

types_alleles = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];

tbl = zeros(6,6,6);

for i = 1:6
    type1 = i;
    for j = 1:6
        type2 = j;
        for k = 1:6
            type3_alleles = types_alleles(k,:);
            if type3_alleles(1) == type3_alleles(2)
                prob = ...
                    gamete_prob(type1,type3_alleles(1),a) * ...
                    gamete_prob(type2,type3_alleles(2),a);
            else
                prob_lr = ...
                    gamete_prob(type1,type3_alleles(1),a) * ...
                    gamete_prob(type2,type3_alleles(2),a);
                prob_rl = ...
                    gamete_prob(type1,type3_alleles(2),a) * ...
                    gamete_prob(type2,type3_alleles(1),a);
                prob = prob_lr + prob_rl; 
            end
            tbl(i,j,k) = prob;
        end
    end
end

end


%=========================================================================%
function prob = gamete_prob(type, gamete, a)

% Returns the probability of a type "type" individual producing a gamete of
% type "gamete", where type is a number between 1 and 6, and a gamete is a
% number between 1 and 3, as described at the top of this file.

P = a.P;
q = a.q;

if gamete == 1
    if type == 1
        prob = 1;
    elseif type == 2
        prob = (1-q)/2;
    elseif type == 3
        prob = 1/2;
    else
        prob = 0;
    end
elseif gamete == 2
    if type == 2
        prob = (1+q*P)/2;
    elseif type == 4
        prob = 1;
    elseif type == 5
        prob = 1/2;
    else
        prob = 0;
    end
elseif gamete == 3
    if type == 2
        prob = q*(1-P)/2;
    elseif type == 3
        prob = 1/2;
    elseif type == 5
        prob = 1/2;
    elseif type == 6
        prob = 1;
    else
        prob = 0;
    end
else
    error('Unknown gamete.')
end

end