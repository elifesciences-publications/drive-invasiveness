function [T, Y, parms] = ...
    simulate_deterministic(strategy, rel_size, n, P, q, r_cost, d_cost, fname, save_bool)

% rel_size    = 0.01;
xDD_init    = rel_size;
x00_init    = 1-rel_size;
xD0_init    = 0;
g           = 0;

if nargin == 0
    strategy    = 'proposed';
    n           = 5;
    d_cost      = 0.45;
    r_cost      = 0.45;
    q           = 0.95;
    P           = 0.95;
    fname       = 'Results';
    save_bool   = 0;
end

%Fitness parameters
f = struct();
if strcmp(strategy,'classic')
    f.aa = 1;
    f.dd = (1-d_cost)*(1-r_cost);
    f.rr = 1-r_cost;
    f.ad = 1-d_cost;
    f.ar = 1;
    f.dr = (1-d_cost)*(1-r_cost);
    
elseif strcmp(strategy,'proposed')
    f.aa = 1;
    f.dd = 1-d_cost;
    f.rr = 1-r_cost;
    f.ad = 1-d_cost;
    f.ar = 1;
    f.dr = 1-d_cost;
    
end

%Simulation parameters
t_max = 200;

%Storage
a_names = gen_allele_lookup_arr(n);                     %allele name array
d_names = gen_genotp_lookup_arr(a_names);               %diploid name array
[parr, parr_names] = gen_p_arr(n, q, P, g, a_names);    %p(IJ,K) array
if sum(abs(sum(reshape(parr,2*n+2,[])) - 1) > 1E-12) > 0
    abs(sum(reshape(parr,2*n+2,[])) - 1) 
    for i = 1:size(parr,1)
        disp([parr_names{i} char(9) num2str(parr(i))])
    end
    error('p probabilities don''t sum to one.')
end
[farr, farr_names] = gen_f_arr(n, f, a_names);  %fitness array and names

%Individual-frequency storage
d_freqs = zeros(length(d_names),1);
init_arr = d_freqs;
init_arr(1) = xDD_init;
init_arr(2) = xD0_init;
init_arr(2*n+3) = x00_init;

% Set up the parms struct to pass to the simulation function
parms = struct();
parms.n = n;
parms.dnames = d_names;
parms.anames = a_names;
parms.farr = farr;
parms.farr_names = farr_names;
parms.parr = parr;
parms.parr_names = parr_names;

%Simulate
% options = odeset('Events',@(t,y) event(t,y,parms), 'RelTol', 1E-8, ...
%     'AbsTol', 1E-6);
options = odeset('RelTol', 1E-8, 'AbsTol', 1E-8);
[T,Y] = ode45(@(t, y) update_function(t, y, parms), [0 t_max], ...
    init_arr, options);

Y = convert_to_allele_frequencies(n, Y);

parms.d_cost = d_cost;
parms.r_cost = r_cost;
parms.P = P;
parms.q = q;
parms.g = g;
parms.t_max = t_max;
parms.x00_init = x00_init;
parms.xD0_init = xD0_init;
parms.xDD_init = xDD_init; 

if save_bool
    save( ...
        [fname '/' strategy ...
        '_n_' num2str(n) ...
        '_r_' num2str(r_cost) ...
        '_d_' num2str(d_cost) ...
        '.mat' ...
        ], ...
        'max_drive_frequency', ...
        'max_drive_frequency_time', ...
        'parms', ...
        'T', ...
        'Y')
end

end

function [value,isterminal,direction] = event(~,Y,parms)
% dY = update_function(t,Y,parms);

% totalDriveFreqChange = dY(1);
totalDriveFreq = Y(1);

% totalDriveFreqChange = totalDriveFreqChange + 0.5 * sum(dY(2:2*parms.n+2));
totalDriveFreq = totalDriveFreq + 0.5 * sum(Y(2:2*parms.n+2));

% value = [totalDriveFreqChange, totalDriveFreq-(5e-3), (1-1e-3)-totalDriveFreq];
% isterminal = [1, 1, 1];
% direction = [-1, -1, -1];

value = [totalDriveFreq-(1e-3), totalDriveFreq-(1-1e-1)];
isterminal = [1, 1];
direction = [-1, -1];

% value = totalDriveFreq-(1E-3);
% isterminal = 1;
% direction = -1;

end

function [farr, farr_names] = gen_f_arr(n, f, alr)

allele_count = 2*n+2;
total_count = (1/2)*allele_count*(allele_count+1);

farr = zeros(total_count,1);
farr_names = cell(total_count,1);

idx = 1;
for i = 1:2*n+2
    for j = i:2*n+2
        str = [alr{i} alr{j}];
        farr_names{idx} = str;

        if strcmp(alr{i}(1),'D')
            if strcmp(alr{j}(1),'D')
                farr(idx) = f.dd;
            elseif strcmp(alr{j}(1),'S')
                farr(idx) = f.ad;
            elseif strcmp(alr{j}(1),'R')
                farr(idx) = f.dr;
            end
        elseif strcmp(alr{i}(1),'S')
            if strcmp(alr{j}(1),'S')
                farr(idx) = f.aa;
            elseif strcmp(alr{j}(1),'R')
                farr(idx) = f.ar;
            end
        elseif strcmp(alr{i}(1),'R')
            if strcmp(alr{j}(1),'R')
                farr(idx) = f.rr;
            end
        end
        idx = idx + 1;
    end
end

end

function [parr, parr_names] = gen_p_arr(n, q, P, g, alr)

allele_count = 2*n+2;
total_count = allele_count*((1/2)*allele_count*(allele_count+1));

parr = zeros(total_count,1);
parr_names = cell(total_count,1);

p_L = @(resist_i, lost_total, cuts, n) (n-resist_i-lost_total+1) * ...
    nchoosek(lost_total-2,cuts-2)/nchoosek(n-resist_i,cuts);
p_C = @(resist_i, cuts, q, n) nchoosek(n-resist_i,cuts)*q^cuts*...
    (1-q)^(n-resist_i-cuts);

idx = 1;
for i = 1:2*n+2
    for j = i:2*n+2
        for k = 1:2*n+2
            str = [alr{i} alr{j} ',' alr{k}];
            parr_names{idx} = str;
            
            %p is only interesting if a drive heterozygote
            if strcmp(alr{i},'D') %D is index 1, so only have to check i
                
                if j == i 
                    %This case represents DD individual
                    parr(idx) = (k == i);
                    
                else
                    %This case is DX where X = S_i or R_i
                    
                    if strcmp(alr{j}(1), 'R')
                        
                        if strcmp(alr{k}(1), 'R')
                            % p(DR, R)
                            
                            if j == k
                                % p(DR_i, R_i)    
                                i_val = str2double(alr{k}(2:end));
                                parr(idx) = 0.5 * (1-q)^(n-i_val);
                                
                            elseif k > j
                                % p(DR_i, R_k)
                                i_val = str2double(alr{j}(2:end));
                                j_val = str2double(alr{k}(2:end));
                                
                                if j_val == i_val+1
                                    parr(idx) = 0.5*(1-P)*p_C(i_val,1,q,n);
                                elseif j_val >= i_val+2
                                    val = 0;
                                    for c = 2:j_val-i_val
                                        val = val + ...
                                            p_L(i_val,j_val-i_val,c,n) *...
                                            p_C(i_val,c,q,n);
                                    end
                                    val = val * (1-P)/2;
                                    parr(idx) = val;
                                end
                            else
                                parr(idx) = 0;
                            end
                            
                        elseif strcmp(alr{k}(1), 'S')
                            % p(DR, S)
                            parr(idx) = 0;
                            
                        elseif strcmp(alr{k}(1), 'D')
                            % p(DR, D)
                            i_val = str2double(alr{j}(2:end));
                            parr(idx) = 0.5 + 0.5*P*(1-(1-q)^(n-i_val));
                            
                        end
                        
                    elseif strcmp(alr{j}(1),'S')
                        
                        if strcmp(alr{k}(1), 'R')
                            % p(DS, R)
                            
                            i_val = str2double(alr{j}(2:end));
                            j_val = str2double(alr{k}(2:end));
                            
                            if j_val == i_val+1
                                parr(idx) = 0.5 * p_C(i_val,1,q,n)*...
                                    (1-P)*(1-g);
                            elseif j_val >= i_val+2
                                val = 0;
                                for c = 2:j_val-i_val
                                    val = val + ...
                                        p_L(i_val,j_val-i_val,c,n) * ...
                                        p_C(i_val,c,q,n);
                                end
                                val = val * (1-P)/2;
                                parr(idx) = val;
                            else
                                parr(idx) = 0; 
                            end
                            
                        elseif strcmp(alr{k}(1), 'S')
                            % p(DS, S)
                            
                            if k == j
                                i_val = str2double(alr{k}(2:end));
                                parr(idx) = 0.5*(1-q)^(n-i_val);
                                
                            elseif k>j
                                i_val = str2double(alr{j}(2:end));
                                j_val = str2double(alr{k}(2:end));
                                if j_val == i_val + 1
                                    parr(idx) = 0.5*p_C(i_val,1,q,n)*...
                                        (1-P)*g;
                                elseif j_val >= i_val + 2
                                    parr(idx) = 0;
                                end
                            else
                                parr(idx) = 0;
                                
                            end
                                
                            
                        elseif strcmp(alr{k}(1), 'D')
                            % p(DS, D)
                            
                            i_val = str2double(alr{j}(2:end));
                            parr(idx) = 0.5 + 0.5*P* ...
                                (1-(1-q)^(n-i_val));     
                        end
                    end
                end
                               
            %Otherwise just standard inheritance    
            else
                if k == j
                    parr(idx) = parr(idx) + 1/2;
                end
                if k == i
                    parr(idx) = parr(idx) + 1/2;
                end
            end
            
            idx = idx + 1;
        end
    end
end

end

function glr = gen_genotp_lookup_arr(larr)  

allele_count = length(larr);
glr = cell((1/2)*allele_count*(allele_count+1),1);
idx = 1;
for i = 1:allele_count
    for j = i:allele_count
        str = [larr{i} larr{j}];
        glr{idx} = str;
        idx = idx+1;
    end
end

end

function alr = gen_allele_lookup_arr(n)

alr = cell(2*n+2,1);
alr{1} = 'D';
for i = 2:2+n
    alr{i} = ['S' sprintf('%d', i-2)];
end
for i = 3+n:2+2*n
    alr{i} = ['R' sprintf('%d', i-2-n)];
end

end

function dot_d_freqs = update_function(~, d_freqs, parms)

% d_freqs(d_freqs<1e-10) = 0;
% d_freqs = d_freqs/sum(d_freqs);

%Load parameters
n = parms.n;
anames = parms.anames;
farr = parms.farr;
parr = parms.parr;

%Make the F_array, F_i stores the F_I where anames(i) is the allele name
F_arr = zeros(size(anames));

%Start by calculating F_D
F_arr(1) = d_freqs(1) * farr(1);

for k = 1:n
    F_arr(1) = F_arr(1) + ...
        parr((2*n+2)*(n+2)+(k-1)*(2*n+2)+1) * ...
        farr(n+k+2) * ...
        d_freqs(n+k+2);
end

for k = 0:n
    F_arr(1) = F_arr(1) + ...
        parr((2*n+2)+(k)*(2*n+2)+1) * ...
        farr(k+2) * ...
        d_freqs(k+2);
end


%Now calculate F_{S_i}
for i = 0:n
    idx = i+2;
    for k = 0:n
        krd = (i == k);
        iprime = min(i,k);
        kprime = max(i,k);
        F_arr(idx) = F_arr(idx) + ...
            (1+krd)/2 * ...
            farr(-(1+iprime)*(iprime-4*(1+n))/2+kprime+1-iprime) * ...
            d_freqs(-(1+iprime)*(iprime-4*(1+n))/2+kprime+1);
    end
    for k = 1:n
        F_arr(idx) = F_arr(idx) + ...
            0.5 * farr(-(1+i)*(i-4*(1+n))/2+n-i+1+k) * ...
            d_freqs(-(1+i)*(i-4*(1+n))/2+n-i+1+k);
    end
    for k = 0:i
        F_arr(idx) = F_arr(idx) + ...
            parr((2*n+2)*(k+1)+i+2) * ...
            farr(k+2) * ...
            d_freqs(k+2);
    end
end

%Now calculate F_{R_i}.
for i = 1:n
    idx = n+2+i;
    for k = 1:n
        iprime = min(i,k);
        kprime = max(i,k);
        krd = (i == k);
        F_arr(idx) = F_arr(idx) + ...
            (1+krd)/2 * ...
            farr(-0.5*(-4+iprime-3*n)*(1+iprime+n)-iprime+kprime+1) * ...
            d_freqs(-0.5*(-4+iprime-3*n)*(1+iprime+n)-iprime+kprime+1);
    end
    for k = 0:n
        F_arr(idx) = F_arr(idx) + ...
            0.5 * ...
            farr(-0.5*(1+k)*(k-4*(1+n))+n+i-k+1) * ...
            d_freqs(-0.5*(1+k)*(k-4*(1+n))+n+i-k+1);
    end
    for k = 1:i
        F_arr(idx) = F_arr(idx) + ...
            parr((2*n+2)*(n+2)+(k-1)*(2*n+2)+2+n+i) * ...
            farr(n+k+2) * ...
            d_freqs(n+k+2);
    end
    for k = 0:i
        F_arr(idx) = F_arr(idx) + ...
            parr((2*n+2)+(k)*(2*n+2)+2+n+i) * ...
            farr(k+2) * ...
            d_freqs(k+2);
    end
end


psi = sum(F_arr);

dot_d_freqs = zeros(size(d_freqs));

idx = 1;
for i = 1:2*n+2
    for j = i:2*n+2
        krd = (i == j);
        dot_d_freqs(idx) = (2-krd) * F_arr(i) * F_arr(j) - psi^2 * d_freqs(idx);
        idx = idx+1;
    end
end

end

function [allele_frequencies] = convert_to_allele_frequencies(n, Y)

freqs = zeros(size(Y,1),2*n+2);

idx = 1;
for i = 1:2*n+2
    for j = i:2*n+2
        freqs(:,i) = freqs(:,i) + (1/2)*Y(:,idx);
        freqs(:,j) = freqs(:,j) + (1/2)*Y(:,idx);
        idx = idx + 1;
    end
end

allele_frequencies = freqs;

end