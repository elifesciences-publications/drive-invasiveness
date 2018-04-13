function [ ] = gen_figure_1E_and_F( )

addpath(genpath('Plotting Utilities'))

% Specify the directory where the data can be found
save_dir = '../Data_Storage/Figure_1E_and_1F/';

% Make a directory to save the data if it doesn't exist
if ~exist(save_dir,'dir')
    mkdir(save_dir) 
end

% Generate the data if it doesn't exist in the directory
P_arr = [0.15, 0.5, 0.9];                   % Homing efficiencies
i_max = 10;                                 % Maximum value on the x-axis
generate_data(save_dir, P_arr, i_max);      % Generate data, if not done

close all; figure('position',[768   324   371   253]);
tight_subplot(1,1,1,[0.15,0.05],[0.11,0.03]); hold on

clr = brewermap(8,'RdBu');
c1 = clr(1,:);
lw = 1;

% Make the distribution plot
[m, md] = make_single_distribution_plot(save_dir, P_arr(2), i_max, 8);
plot(1:i_max, md, 's-', 'color', c1, 'linewidth', lw, 'markerfacecolor', c1)
set(gca,'tickdir','out')
set(gca,'ticklength',[0.015, 0.025])
set(gca,'yminortick','on')
set(gca,'xtick',1:i_max)

% Now make the line plots
make_line_plots(save_dir, P_arr, i_max, clr(3:-1:1,:), lw);

end


function [] = generate_data(save_dir, P_arr, i_max)

idx = 1;
for p = 1:length(P_arr)
    for i = 1:i_max
        file_name = [save_dir 'sim_data_P_' num2str(P_arr(p)) ...
            '_i_' num2str(i) '.mat'];
        if ~exist(file_name,'file')
            [a, store_dist] = moran_wm(i, P_arr(p), 1000, inf, 'max_dist');
            save(file_name, 'a', 'store_dist')
        end
        disp(['File ' num2str(idx) ' out of ' ...
            num2str(length(P_arr)*i_max) ' exists or was generated.'])
        idx = idx + 1;
    end
end

end


function [ ] = make_line_plots(save_dir, P_arr, i_max, clrs, lw)

arr_mean = zeros(i_max,length(P_arr));
arr_med = zeros(i_max,length(P_arr));
for p = 1:length(P_arr)
    for i = 1:i_max
        f=load([save_dir 'sim_data_P_' num2str(P_arr(p)) '_i_' num2str(i) '.mat']);
        dist = f.store_dist;
        arr_mean(i,p) = mean(dist);
        arr_med(i,p) = median(dist);
    end
end

figure('position',[ 768   324   371   253]);
tight_subplot(1,1,1,[0.15,0.05],[0.11,0.03]); hold on
for p = 1:length(P_arr)
    plot(1:i_max, arr_mean(:,p), '.-', 'color', clrs(p,:), 'linewidth', lw, 'markersize', 18)
    plot(1:i_max, arr_med(:,p), 's:', 'color', clrs(p,:), 'linewidth', lw, 'markerfacecolor', clrs(p,:))
end

xlabel('Release size')
ylabel('Maximum drive frequency')
set(gca,'box','off')
set(gca,'ticklength',[0.015, 0.025])
set(gca,'yminortick','on')
set(gca,'tickdir','out')
set(gcf,'color','w')

end


function [arr_mean, arr_med] = make_single_distribution_plot(save_dir, P, plots, clr_idx)

dist_cell = cell(plots,1);
arr_mean = zeros(plots,1);
arr_med = zeros(plots,1);

for i = 1:plots
    f=load([save_dir 'sim_data_P_' num2str(P) '_i_' num2str(i) '.mat']);
    dist = f.store_dist;
    dist_cell{i} = dist;
    arr_mean(i) = mean(dist);
    arr_med(i) = median(dist);
end

clr = brewermap(20,'RdBu');

violin(dist_cell'               , ...
    'edgecolor','none'          ,...
    'facecolor',clr(clr_idx,:)  , ...
    'facealpha', 1              , ...
    'mc', [], ...
    'medc', [], 'bw', 0.02);

xlabel('Release size')
ylabel('Maximum drive frequency')
set(gca,'box','off')
set(gca,'tickdir','out')
set(gca,'ticklength',[0.1, 0.1])
set(gcf,'color','w')

end