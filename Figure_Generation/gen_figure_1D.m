function [ ] = gen_figure_1D( )

addpath(genpath('Plotting Utilities'))
addpath(genpath('Simulation Functions'))

% Run 1000 simulations for each of three different scenarios:
%   P = 0.15
%   P = 0.5
%   P = 0.9
%
% Then plot the allele frequencies, showing deterministic results,
% means of simulations, and 50% confidence regions.
% 
% Assuming release size of 15 for each scenario.

% In each of the below, "a" is a struct holding the parameters, and "c" is
% a cell array storing full trajectories for each of the simulations
[a1, c1] = moran_wm(15, 0.15, 1000, 80, 'full');
[a2, c2] = moran_wm(15, 0.50, 1000, 80, 'full');
[a3, c3] = moran_wm(15, 0.90, 1000, 80, 'full');

% Make the figure
figure('position',[509   336   695   224], 'color', 'w');
axs = tight_subplot(1,3,0.05,[0.18,0.05],[0.06,0.02]);

% Plot each of the panels separately
plot_simulations_dist(a1,c1,axs(1));
plot_simulations_dist(a2,c2,axs(2)); set(gca,'ylabel',[])
plot_simulations_dist(a3,c3,axs(3)); set(gca,'ylabel',[])

end


function [ ] = plot_simulations_dist(a, store_cell, ax)

clr         = brewermap(6,'RdBu');
sims        = a.sims;
gens        = a.max_steps / a.N;
sim_alpha   = 1;
face_alpha  = 1;
plot_example_sims   = 1;

% Deterministic simulation
[T, Y] = simulate_deterministic('classic',a.i/a.N,1,a.P,a.q,a.s,a.c,[],0);
D_deterministic = Y(:,1);
R_deterministic = Y(:,3)+Y(:,4);

% Find the longest simulation
max_len = a.max_steps+1;
xvals = (0:a.max_steps)/a.N;

% Extract abundances of D's and R's over time and store
sim_plot_arr_D = zeros(max_len, plot_example_sims);
sim_plot_arr_R = zeros(max_len, plot_example_sims);
save_plot_idxs = randperm(sims); save_plot_idxs = save_plot_idxs(1:plot_example_sims);

saved_total = 0;
all_D = zeros(max_len, sims);
all_R = zeros(max_len, sims);
for i = 1:sims
    
    % Get simulations from storage
    vec = store_cell{i};
    D = vec(:,1);
    R = vec(:,2);
    
    if length(D) < max_len
        all_D(1:length(D),i) = D;
        all_D(length(D)+1:end,i) = D(end)*ones(max_len-length(D),1);
    else
        all_D(:,i) = D(1:max_len);
    end
    if length(R) < max_len
        all_R(1:length(R),i) = R;
        all_R(length(R)+1:end,i) = R(end)*ones(max_len-length(R),1);
    else
        all_R(:,i) = R(1:max_len);
    end
    
    % Save some example plots later
    if sum(i == save_plot_idxs) > 0
        sim_plot_arr_D(:,saved_total+1) = all_D(:,i);
        sim_plot_arr_R(:,saved_total+1) = all_R(:,i);
        saved_total = saved_total+1;
    end
end

% Normalize
all_D = all_D / (2*a.N);
all_R = all_R / (2*a.N);
sim_plot_arr_D = sim_plot_arr_D / (2*a.N);
sim_plot_arr_R = sim_plot_arr_R / (2*a.N);

% Get the means
avg_drive = mean(all_D,2);
avg_resist = mean(all_R,2);

% Get the 95% bounds
bounds_D = quantile(all_D,[0.25 0.75], 2);
bounds_R = quantile(all_R,[0.25 0.75], 2);
interp_x = 0:.1:max(xvals);
interp_bounds_D = interp1(xvals,bounds_D,interp_x);
interp_bounds_R = interp1(xvals,bounds_R,interp_x);

% Make the figure
axes(ax); hold on;

% R -- Make the patches
fl_R = patch([interp_x, interp_x(end:-1:1)], [interp_bounds_R(:,1)', ...
    interp_bounds_R(end:-1:1,2)'], clr(4,:));
set(fl_R,'facecolor',clr(4,:))
set(fl_R,'linestyle','none')
set(fl_R,'facealpha',face_alpha)

% D -- Make the patches
fl_D = patch([interp_x, interp_x(end:-1:1)], [interp_bounds_D(:,1)', ...
    interp_bounds_D(end:-1:1,2)'], clr(3,:));
set(fl_D,'facecolor',clr(3,:))
set(fl_D,'linestyle','none')
set(fl_D,'facealpha',face_alpha)

if plot_example_sims > 0
    % R -- Individual simulations
    for i = 1:plot_example_sims
        p2 = plot(xvals,sim_plot_arr_R(:,i),'Color',clr(4,:),'linewidth',1);
        p2.Color(4) = sim_alpha;
    end
end

% R -- Deterministic
plot(ax,T,R_deterministic,'-','color',clr(6,:),'linewidth',1.5);

% R -- Average
plot(ax,xvals,avg_resist,'-','color',clr(5,:),'linewidth',1.5);

if plot_example_sims > 0
    % D -- Individual simulations
    for i = 1:plot_example_sims
        p1 = plot(ax,xvals,sim_plot_arr_D(:,i),'Color',clr(3,:),'linewidth',1);
        p1.Color(4) = sim_alpha;
    end
end

% D -- Deterministic
plot(ax,T,D_deterministic,'-','color',clr(1,:),'linewidth',1.5);

% D -- Plot the averages
plot(xvals,avg_drive,'-','color',clr(2,:),'linewidth',1.5);

% Format the plot
set(ax,'TickDir','Out')
xlim([0 gens])
xlabel('Time (generations)')
ylabel('Allele frequency')
set(gca,'ticklength',[0.025 0.025])
set(gca,'XMinorTick','on','YMinorTick','on')

end