function [ ] = gen_figure_2E( )

addpath(genpath('Plotting Utilities/'))
addpath('Simulation Functions/')

cm = brewermap(4,'Reds');

data_directory = '../Data_Storage/Figure_2E/';

n = 10000; % simulations
m_arr = logspace(-5, log10(0.5), 51)';
esc_pop_arr = (1:4)';

prob_arr = zeros(4,length(m_arr));
for i = 0:50
    try
        vec = load([data_directory 'data_index_' num2str(i) '.csv']);
        prob_arr(:,i+1) = vec;
    catch
        prob_arr(:,i+1) = zeros(4,1);
    end
end

close all
figure('position',[856   413   306   222]);
leg_arr = zeros(4,1);
for i = 1:length(esc_pop_arr)
    if i == 1; hold on; end
    leg_arr(i) = plot(m_arr, prob_arr(i,:), '.-','color', cm(i,:));
    
    for j = 1:length(m_arr)
        val = prob_arr(i,j);
        st_er = sqrt(val*(1-val)/n);
        plot([m_arr(j), m_arr(j)], [val-st_er, val+st_er], ...
            'color', cm(i,:));
    end
    
end
set(gca,'xscale','log')
set(gca,'position',[ 0.1088    0.1628    0.6484    0.7622])
set(gca,'tickdir','out')
xlim([min(m_arr), max(m_arr)])
axis square
set(gca,'xtick',[1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
xlabel('Migration rate')
set(gca,'ticklength',[.02,.02])
ylabel('Escape probability')
legend(leg_arr,{'1 population','2 populations','3 populations','4 populations'})

end

%=========================================================================%
function [plot_cell, count_cell] = gather_data(fpath, m_arr, esc_pop_arr)

plot_cell = cell(length(esc_pop_arr),1);
count_cell = cell(length(esc_pop_arr),1);
for i = 1:length(plot_cell)
    plot_cell{i} = zeros(length(m_arr),1); 
    count_cell{i} = zeros(length(m_arr),1);
end

d = dir(fpath);
d = d(3:end);

for i = 1:length(d)
    f = load([fpath 'data_' num2str(i) '.mat']);
    try
        s = f.sim;
    catch
        s = f.sims;
    end
    escape_prob = sum(f.escape_arr(1:s)) / s;
    [m, esc_pop] = get_parameters(i, m_arr, esc_pop_arr);
    idx_m = find(m_arr == m);
    idx_ep = find(esc_pop_arr == esc_pop);
    plot_cell{idx_ep}(idx_m) = escape_prob;
    count_cell{idx_ep}(idx_m) = s;
end


end


%=========================================================================%
function [m, esc_pop] = get_parameters(idx, m_arr, esc_pop_arr)

len_m = length(m_arr);
len_esc_pop = length(esc_pop_arr);

param_arr = zeros(len_m * len_esc_pop, 2);

if idx > len_m * len_esc_pop
    error('Requested parameter index invalid.') 
end

param_arr(:,2) = repmat(esc_pop_arr,len_m,1);
for i = 1:len_m
    j = len_esc_pop * (i - 1) + 1;
    param_arr(j:j+len_esc_pop-1,1) = m_arr(i);
end

m = param_arr(idx, 1);
esc_pop = param_arr(idx, 2);

end