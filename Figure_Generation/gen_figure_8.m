function [ ] = gen_figure_8( )

addpath(genpath('Plotting Utilities/'))

fname = '../Data_Storage/Figure_8/data.csv';

close all
figure('units','normalized','position',[0.6141    0.3975    0.2805    0.3913]);
tight_subplot(1,1,1,[0.15,0.05],[0.11,0.03]);

cmap = brewermap(100,'*rdbu');
arr = csvread(fname);
arr(arr<0) = nan;
arr(:,end) = [];

c_vals = 0:0.01:0.5;
f_vals = 1-c_vals;

s_vals = 0:0.02:1;
res_fit_vals = 1-s_vals;

imagesc(res_fit_vals,f_vals,arr')
set(gca,'clim',[15/500,1])
axis square
colormap(cmap)
set(gca,'ydir','normal')
colorbar
xlabel('Resistance fitness, 1-s')
ylabel('Drive fitness, f')
set(gca,'tickdir','out')
set(gca,'box','off')
set(gca,'ytick',0.5:.1:1)
set(gca,'yminortick','on')
set(gca,'xminortick','on')
set(gca,'ticklength',[0.018, 0.018])

end