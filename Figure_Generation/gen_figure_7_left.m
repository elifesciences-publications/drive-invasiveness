function [ ] = gen_figure_7_left( )

addpath(genpath('Plotting Utilities/'))

fname = '../Data_Storage/Figure_7_left/data.csv';

close all
figure('units','normalized','position',[0.6141    0.3975    0.2805    0.3913]);
tight_subplot(1,1,1,[0.15,0.05],[0.11,0.03]);

cmap = brewermap(100,'*rdbu');
arr = csvread(fname);
arr(arr<0) = nan;
arr(:,end) = [];

c_vals = 0:0.01:0.5;
f_vals = 1-c_vals;

P_vals = 0:.02:1;
imagesc(P_vals,f_vals,arr')
set(gca,'clim',[15/500,1])
axis square
colormap(cmap)
set(gca,'ydir','normal')
colorbar
xlabel('Homing efficiency, P')
ylabel('Drive fitness, f')
set(gca,'tickdir','out')
set(gca,'box','off')
set(gca,'ytick',0.5:.1:1)
set(gca,'yminortick','on')
set(gca,'xminortick','on')
set(gca,'ticklength',[0.018, 0.018])

hold on;
plot(P_vals,1./(1+P_vals),'w','linewidth',1)

effs = 0.01*[
    100
62
52
98
14
83
95
99
];
ylim([0.5 1])

hold on
for i = 1:length(effs)
    e = effs(i);
    plot([e, e], [1, 1/(1+e)], ':', 'color', 'w', 'linewidth', 1.5);
end

end