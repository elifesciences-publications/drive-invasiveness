function [ ] = gen_figure_2D( )

addpath(genpath('Plotting Utilities/'))
addpath('Simulation Functions/')

data_directory = '../Data_Storage/Figure_2D/';

close all
prob_arr = zeros(51,51);
for i = 0:50
    try
        vec = load([data_directory 'data_index_' num2str(i) '.csv']);
        prob_arr(i+1,:) = vec;
    catch
        prob_arr(i+1,:) = zeros(1,51);
    end
end

m_arr = logspace(-5,log10(0.5),size(prob_arr,2));
P_arr = linspace(0,1,size(prob_arr,1));

cm = brewermap(100,'Reds');

% f = gcf;
f=figure('position',[856   413   306   222]);
set(f,'color','w')

plt=surf(m_arr, P_arr, prob_arr);
set(plt,'linestyle','none')
colormap(cm)

set(gca,'xscale','log')
set(gca,'xlim',[min(m_arr) max(m_arr)])
view(0,90)
xlabel('Migration rate')
ylabel('Homing efficiency')
colorbar
axis square
set(gca,'clim',[0 1])

set(gca,'xtick',[1e-5, 1e-4, 1e-3, 1e-2, 1e-1])

set(gca,'tickdir','out')
set(gca,'box','on')
set(gca,'ticklength',[.02,.02])

end