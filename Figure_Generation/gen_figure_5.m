function [ ] = gen_figure_5( )

addpath(genpath('Plotting Utilities'))
folder = '../Data_Storage/Figure_5/';
fnames = dir(folder);

% Names are of the form: MaxFreqDist_r_0.04.csv
bs_arr = get_variable_vals(fnames);
dist_cell = cell(1,length(bs_arr));
for bs_i = 1:length(bs_arr)
    dist = csvread([folder '/data_k_' num2str(bs_arr(bs_i)) '.csv']);
    dist_cell{bs_i} = dist;
end

addpath('..')
close all; figure('units','normalized','position',[0.5000    0.2775    0.3570    0.5112]);
ha=tight_subplot(2,1,0.08,[0.1,0.03],[0.1,0.03]);
h_diff = 0.07;
pos1=get(ha(1),'position');
pos2=get(ha(2),'position');
pos1=pos1+[0,h_diff,0,-h_diff];
pos2=pos2+[0,0,0,h_diff];
set(ha(1),'position',pos1)
set(ha(2),'position',pos2)

clrs = brewermap(6,'RdBu');
lw = 1;

axes(ha(2)); hold on;
violin(dist_cell               , ...
    'edgecolor','none'          ,...
    'facecolor',clrs(3,:)  , ...
    'facealpha', 1              , ...
    'mc', [], ... %clr(1,:)              , ...
    'medc', [], 'bw', 0.02);
hold on;
yl = get(gca,'ylim');
yl(2) = 0.65;
set(gca,'ylim',yl)

f = @(x) median(x);
p2 = plot(1:length(bs_arr), cellfun(@mean,dist_cell), '.-', 'color', ...
    clrs(1,:), 'linewidth', lw, 'markersize', 18);
% p3 = plot(1:length(bs_arr), cellfun(f,dist_cell), '-x', 'color', ...
%     clrs(1,:),'linewidth', lw, 'markersize', 6);
set(gca,'xtick',1:length(bs_arr))
set(gca,'xticklabel',strread(num2str(bs_arr),'%s'));
ylabel('Maximum drive frequency')
xlabel('Number of offspring per mating')
set(gca,'tickdir','out')
ch = get(gca,'children');
set(gca,'box','off')
p1 = ch(end);
xl = get(gca,'xlim');
leg = legend([p1,p2],{'Distribution','Mean'});
set(leg,'edgecolor','none')
set(leg,'location','southwest')

axes(ha(1));
k = 1:25;
plot(k,250*(2*k+6)/4,'.-','color',clrs(5,:),'markersize',18); hold on;
plot(k,8*(2*k+6)/4,'.-','color',clrs(6,:),'markersize',18);
set(gca,'yscale','log')
set(gca,'xlim',xl)
set(gca,'xtick',1:length(bs_arr))
set(gca,'xticklabel',strread(num2str(bs_arr),'%s'))
set(gca,'tickdir','out')
set(gca,'box','off')
ylabel('Individuals')
leg=legend('Total population','Release size');
set(leg,'location','northwest')
set(leg,'edgecolor','none')
ylim([10,10000])
ax = gca;
ax.YGrid = 'on';

set(gcf,'color','w')

end

function r_arr = get_variable_vals(fnames)

for i = length(fnames):-1:1
    fname = fnames(i).name;
    if fname(1) == '.'
        fnames(i) = [];
    elseif isdir([fnames(i).folder '/' fname])
        fnames(i) = []; 
    end
end

r_arr = zeros(1,length(fnames));
for i = 1:length(fnames)
    f = fnames(i).name;
    fs = strsplit(f,{'_'});
    r_arr(i) = str2double(fs{3}(1:end-4));
end

r_arr = unique(r_arr);

end