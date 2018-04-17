function [ ] = gen_figure_1B( )

addpath(genpath('Plotting Utilities'))

% Load colors from brewermap
clr = brewermap(3,'Reds');
clr1 = clr(1,:);
clr2 = clr(end,:);

% Make the figure
close all
figure('position',[817   282   256   241],'color','w')
p_vals = 0:.01:1;
f_vals = 1 ./ (1+p_vals);
p=patch([p_vals,1,0],[f_vals,1,1],clr1);
uistack(p,'bottom')
set(p,'edgecolor','none')
xlabel('Homing efficiency')
ylabel('Heterozygote fitness')
set(gca,'tickdir','out')
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'ticklength',[0.025 0.025])
set(gca,'layer','top')
set(gca,'box','on')

% Plot the homing efficiencies from Table S1.
effs = 0.01*[100, 62, 52, 98, 14, 83, 95, 99];
ylim([0.5 1])

hold on
for i = 1:length(effs)
    e = effs(i);
    plot([e, e], [1/(1+e), 1], 'color', clr2, 'linewidth', 1);
end

end