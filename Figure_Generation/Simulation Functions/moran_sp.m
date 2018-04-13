function [varargout] = moran_sp(m, subpops, P, i_0, N, sims, quitearly)

close all

% Simulates the Moran process with sexual reproduction and three alleles: 
% WT, drive, and resistant. This function considers subpopulations.

% In vectors used throughout, type ordering:
%   1 - WW
%   2 - WD
%   3 - WR
%   4 - DD
%   5 - DR
%   6 - RR
%
% Gamete ordering:
%   1 - W
%   2 - D
%   3 - R

keep_sims = 1;

a.q = 1.0;  % cutting
a.c = 0.1;  % fitness drive
a.s = 0.0;  % fitness resistance

if nargin == 0
    a.m = 2^-6;     % migration rate
    sims = 500;     % Number of simulations
    a.S = 5;        % # of subpopulations
    a.P = 0.5;      % homing efficiency
    a.i = 10;       % initial DDs
    a.N = 100;  % subpopulation size
else
    a.S = subpops;
    a.m = m;
    a.P = P;
    a.i = i_0;
    a.N = N;
end

G = mig_arr(a);             % migration graph
tbl = mating_table(a);      % mating table
birth = init_birth(a);      % birth rates
death = init_death(a);      % death rates

if keep_sims
    store_cell = cell(sims,1); 
end

max_arr = zeros(sims,1);
for sim = 1:sims
    pop = init_pop(a);
    
    if keep_sims
        store_vec = zeros(6,a.S,1000000);
        store_vec(:,:,1) = pop;
        len = 1;
    end
    
    max_d_freq = 0;
    while sum(sum(pop,2) == a.N) == 0
        pop = update(pop,tbl,birth,death,a,G);
        cts = allele_counts(pop);
        d_ct = sum(cts(2,:),2) / (2*a.N);
        if d_ct > max_d_freq
            max_d_freq = d_ct; 
        end
        if keep_sims
            store_vec(:,:,len+1) = pop; 
            len = len + 1;
        end
        
        if quitearly && (sum(pop(4,:) + pop(2,:,:) + pop(5,:,:)) == 0)
            break 
        end
    end
    if keep_sims
        store_vec(:,:,len+1:end) = [];
        store_cell{sim} = allele_counts(store_vec);
    end
    max_arr(sim) = max_d_freq;
    
    % Message and save occassionally
    if ~mod(sim,10)
        disp(['Done with ' num2str(sim) ' sims']);
    end
end

if keep_sims
    varargout{1} = store_cell;
end

end


%=========================================================================%
function cts = allele_counts(pop)

cts = zeros(3,size(pop,2),size(pop,3));
cts(1,:,:) = 2*pop(1,:,:) + pop(2,:,:) + pop(3,:,:);
cts(2,:,:) = 2*pop(4,:,:) + pop(2,:,:) + pop(5,:,:);
cts(3,:,:) = 2*pop(6,:,:) + pop(3,:,:) + pop(5,:,:);

end


%=========================================================================%
function arr = mig_arr(a)

% Complete graph
arr = ones(a.S);
% arr = gallery('tridiag',ones(1,a.S-1),ones(1,a.S),ones(1,a.S-1));
% arr = full(arr);

end


%=========================================================================%
function pop_out = update(pop_in, tbl, birth, death, a, G)

if rand <= a.m
    % Send off to the migration function
    pop_out = migrate(pop_in, a, G);
else
    % Pick a subpopulation proportional to total fitness
    total_fitness = sum(pop_in .* repmat(birth,1,a.S),1);
    sp_idx = find(rand <= cumsum(total_fitness/sum(total_fitness)),1);
    
    % Update it
    subpop = update_single_pop(pop_in(:,sp_idx), tbl, birth, death);
    pop_out = pop_in;
    pop_out(:,sp_idx) = subpop;
end

end


%=========================================================================%
function pop_out = migrate(pop_in, a, G)

% Pick a subpopulation proportional to size
subpop1 = find(rand <= cumsum(sum(pop_in,1)/a.N), 1);

% Pick an individual uniformly from the subpopulation
indiv = find(rand <= cumsum(pop_in(:,subpop1)/sum(pop_in(:,subpop1))),1);

% Pick the destination subpopulation based on graph
idxs = find(G(subpop1,:)>0);
idxs(idxs == subpop1) = [];
subpop2 = idxs(find(rand <= cumsum(G(subpop1,idxs)/sum(G(subpop1,idxs))), 1));

% Move the individual
pop_out = pop_in;
pop_out(indiv,subpop1) = pop_out(indiv,subpop1)-1;
pop_out(indiv,subpop2) = pop_out(indiv,subpop2)+1;

end


%=========================================================================%
function pop_out = update_single_pop(pop_in, tbl, birth, death)

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

pop = zeros(6,a.S);
ww_each = a.N / a.S;
pop(1,2:end) = ww_each;
pop(1,1) = ww_each - a.i;
pop(4,1) = a.i;

end


%=========================================================================%
function fit = init_death(~)

fit = ones(6,1);

end


%=========================================================================%
function fit = init_birth(a)

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