function [ ] = gen_figure_9( )

addpath(genpath('Plotting Utilities/'))

folder = '../Data_Storage/Figure_9/';
fnames = dir(folder);

[F_arr, P_arr] = get_variable_vals(fnames);
dist_cell = cell(length(F_arr),length(P_arr));
for F_i = 1:length(F_arr)
    for P_i = 1:length(P_arr)
        fname = [folder '/data_P_' num2str(P_arr(P_i)) '_F_' num2str(F_arr(F_i)) '.csv'];
        if exist(fname,'file')
            dist = csvread([folder '/data_P_' num2str(P_arr(P_i)) '_F_' num2str(F_arr(F_i)) '.csv']);
        else
            dist = 0;
        end
        dist_cell{F_i,P_i} = dist;
    end
end

close all; figure('units','normalized','position',[0.5703    0.2950    0.3469    0.4938]);
ha = tight_subplot(3,1,.065,[0.1,0.03],[0.09,0.025]);
set(gcf,'color','w')

clrs = brewermap(6,'RdBu');

yl_min = inf;
yl_max = -inf;

for i = 1:3
    axes(ha(3-i+1))
    violin(dist_cell(:,i)'               , ...
        'edgecolor','none'          ,...
        'facecolor',clrs(3,:)  , ...
        'facealpha', 1              , ...
        'mc', [], ... %clr(1,:)              , ...
        'medc', [], 'bw', 0.02);
        hold on;
    set(gca,'xtick',1:5:51)
    set(gca,'xticklabel',strread(num2str(0:.1:1),'%s'))
    
    plot(cellfun(@mean,dist_cell(:,i)),'.-','markersize',6,'color',clrs(1,:),'linewidth',0.5)
    yl = get(gca,'ylim');
    if yl(1) < yl_min
        yl_min = yl(1);
    end
    if yl(2) > yl_max
        yl_max = yl(2); 
    end
    if i == 1
        xlabel('Selfing rate, s') 
    end
    if i == 3
        ch = get(gca,'children');
        l1 = ch(1);
        l2 = ch(2);
        leg = legend([l1, l2],{'Mean','Distribution'});
        set(leg,'location','southwest')
        set(leg,'edgecolor','none')
    end
    ylabel('Max drive freq')
    set(gca,'tickdir','out')
    set(gca,'box','off')
    text(45.5,0.9,['P = ' sprintf('%1.2f',P_arr(i))])
end
for i = 1:3
    set(ha(i),'ylim',[yl_min yl_max]) 
end

end

function [F_arr, P_arr] = get_variable_vals(fnames)

for i = length(fnames):-1:1
    fname = fnames(i).name;
    if fname(1) == '.'
        fnames(i) = [];
    elseif isdir([fnames(i).folder '/' fname])
        fnames(i) = [];
    end
end

P_arr = zeros(1,length(fnames));
F_arr = zeros(1,length(fnames));
for i = 1:length(fnames)
    f = fnames(i).name;
    fs = strsplit(f,{'_'});
    P_arr(i) = str2double(fs{3});
    F_arr(i) = str2double(fs{5}(1:end-4));
end

F_arr = unique(F_arr);
P_arr = unique(P_arr);

end