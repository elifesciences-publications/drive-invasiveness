function [ ] = gen_figure_2C( )

addpath(genpath('Plotting Utilities/'))
addpath('Simulation Functions/')

m_arr = [10^-1, 10^-2, 10^-3, 10^-4];
subpops = 5;
P       = 0.5;
i_0     = 15;
N       = 500;
sims    = 1;

store_cell_cell = cell(length(m_arr),1);
max_len = 100;
for i = 1:length(m_arr)
    store_cell_cell{i} = moran_sp(m_arr(i), subpops, P, i_0, N, sims, 1);
    len = find(sum(squeeze(store_cell_cell{i}{1}(2,:,:)),1)==0,1,'first');
    if len > max_len
        max_len = len; 
    end
end

close all
clrs = brewermap(subpops,'Reds');
figure('position',[855   319   349   262]);
set(gcf,'color','w')
ha = tight_subplot(length(m_arr),1,0.09,[0.15,0.05],[0.12,0.03]);

for i = 1:length(m_arr)
    single_plot(ha(i), clrs, store_cell_cell{i}, N, max_len)
end

yl_max = 0;
for i = 1:length(m_arr)
    yl = get(ha(i),'ylim');
    if yl(2) > yl_max
        yl_max = yl(2); 
    end
end
ylims = [0, yl_max];
for i = 1:length(m_arr); set(ha(i),'ylim',ylims); end;

axes(ha(end))
xlabel('Time (generations)')

end

function [ ] = single_plot(ax, clrs, store_cell, N, max_len)

vec = store_cell{1};

axes(ax); hold on;

len = size(vec,3);
xvals = (0:len-1) / N;
xinterp = 0:0.1:max(xvals);

% rank the simulations in order of their timing
rank_arr = zeros(size(vec,2),1);
for i = 1:size(vec,2)
    D = squeeze(vec(2,i,:));
    max_D = max(D);
    rank_arr(i) = find(D >= 0.5*max_D,1,'first');
end
[~, idxs] = sort(rank_arr,'ascend');

for i = 1:size(vec,2)
    D = squeeze(vec(2,i,:));
    total = squeeze(sum(vec(:,i,:),1));
    D_interp = interp1(xvals, D./total, xinterp);
    plot(xinterp, D_interp, 'color', clrs(idxs(i),:));
end

xlim([0, (max_len-1) / N])
ylabel('Drive freq.')
set(gca, ...
    'tickdir', 'out', ...
    'ticklength', [0.015, 0.015], ...
    'yminortick', 'on', ...
    'xminortick', 'on')

end