%% Analysing WGCNA output using custom w (wgcna) class (Kevin Coffey & Russell Marx - 2019)
% You must be in the main code directory to run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ALWAYS BACKUP YOUR WGCNA OUTPUT FILES BEFORE RUNNING! %%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run IP Samples
clear all

%% Initialize class and load WGCNA Output
w = WGCNA;

%% Run on first pass to save better colors and smaller 'dissTOM.mat' 
% w.loadGeneTable("\\10.159.50.8\kcoffey\Coffey.Lab\Content.Production\Collaborator.Content\2025.04.01 Bergkamp WGCNA Fresh\WGCNA Output\Gene Modules IP Raw.csv",'coolColors', 1);
% w.loadDissTOM("\\10.159.50.8\kcoffey\Coffey.Lab\Content.Production\Collaborator.Content\2025.04.01 Bergkamp WGCNA Fresh\WGCNA Output\dissTOM IP.csv",'saveAsMat', 1);
% w.colorMap='crayola';

%% Run on all subsequent  
w.loadGeneTable(".\DataSheets\Gene Modules IP.csv",'coolColors', 0);
w.loadDissTOM(".\DataSheets\dissTOM IP.mat",'saveAsMat', 0);
w.colorMap='crayola';

%% Figure and Results Directory
cd('.\Raw Figures');

%% Plot the eigen-genes
% w.createGraphTOM; Only for Weighted 
w.calculateEigenGenes
w.plotEigenGenes
savefig(gcf,'eigenGenes_IP_Raw.fig');
export_fig('eigenGenes_IP_Raw.png','-m3');

%% Merge Eigen-genes
% Don't save colors until you are happy with merging

% Don't run at all if the modules are 
% % w.mergeEigenGenes('mergeThreshold', 1,'saveTable',1);
% % w.calculateEigenGenes
% % w.plotEigenGenes
% % export_fig('eigenGenes_merged.png','-m3')

%% Module Significance w/ DESEQ2
w.loadDESEQ("..\DESeq2\1 DESeq2 Air0 IP vs Air8 IP.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Air Module Significance.png','-m3')

w.loadDESEQ("..\DESeq2\2 DESeq2 EtOH0 IP vs All Air IP.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
savefig(gcf,'EtOH Module Significance.fig');
export_fig('EtOH Module Significance_Tall.png','-m3')

w.loadDESEQ("..\DESeq2\3 DESeq2 ETOH8 IP vs ETOH0 IP.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
savefig(gcf,'Withdrawal Module Significance.fig');
export_fig('Withdrawal Module Significance_Tall.png','-m3')

w.loadDESEQ("..\DESeq2\4 DESeq2 All IP vs All IN.xlsx");
w.plotWaldStatForEachModule('Wald_Stats')
export_fig('Enr Module Significance.png','-m3')

%% MinSpanTree - Module Colors & ETOH Significance
% Create Graph
w.loadDESEQ("..\DESeq2\2 DESeq2 EtOH0 IP vs All Air IP.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM,'Type','forest');
% g = w.removeDisconnectedNodes(g);
% Figure
figure('Color','k','Position',[1 1 450 450]);
h = plot(g)
layout(h,'force','Iterations',120,'WeightEffect','inverse','UseGravity','on');
%layout(h,'layered');
[colorNames, rgb] = colornames(w.colorMap,categories(removecats(g.Nodes.module)));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
h.EdgeColor = [0.9,0.6,0.10];
h.LineWidth = 1.5;
savefig(gcf,'MinSpan Module Color_IP.fig');
export_fig('MinSpan Module Color_IP.png','-m3');
% Add ETOH Significance
h.NodeCData = g.Nodes.Wald_Stats;
set(gcf,'Colormap',crameri('Vanimo'),'Color','k');
c=colorbar('Color','w');
c.Label.String="Wald-Statitic";
c.Label.FontSize=14;
c.FontSize=14;
caxis([-8 8]);
set(gcf,'Position',[1000 100 520 450]);
savefig(gcf,'MinSpan Module Color_ETOH.fig');
export_fig('MinSpan Module Color_ETOH.png','-m3');


%% MinSpanTree - Module Colors & Withdrawal Significance
% Create Graph
w.loadDESEQ("..\DESeq2\3 DESeq2 ETOH8 IP vs ETOH0 IP.xlsx");
w.createGraphTOM;
g = minspantree(w.graphTOM,'Type','forest');
% g = w.removeDisconnectedNodes(g);
% Figure
figure('Color','k','Position',[1 1 450 450]);
h = plot(g)
layout(h,'force','Iterations',120,'WeightEffect','inverse','UseGravity','on');
%layout(h,'layered');
[colorNames, rgb] = colornames(w.colorMap,categories(removecats(g.Nodes.module)));
set(gcf,'Colormap',rgb);
h.NodeCData = findgroups(g.Nodes.module);
axis off
box off
h.EdgeColor = [0.9,0.6,0.10];
h.LineWidth = 1.5;
% Add Withdrawal Significance
h.NodeCData = g.Nodes.Wald_Stats;
set(gcf,'Colormap',crameri('Vanimo'),'Color','k');
c=colorbar('Color','w');
c.Label.String="Wald-Statitic";
c.Label.FontSize=14;
c.FontSize=14;
caxis([-6 6]);
set(gcf,'Position',[1000 100 520 450]);
savefig(gcf,'MinSpan Module Color_Withdrawal.fig');
export_fig('MinSpan Module Color_Withdrawal.png','-m3');

% % 3D Layout for Dope Rotation
% figure('Color','k','Position',[1 1 850 900]);
% h = plot(g,'MarkerSize',3,'EdgeAlpha',.25);
% layout(h,'force3','Iterations',30,'WeightEffect','inverse','UseGravity','off');
% [colorNames, rgb] = colornames(w.colorMap,categories(g.Nodes.module));
% set(gcf,'Colormap',rgb);
% h.NodeCData = findgroups(g.Nodes.module);
% axis off
% box off
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=30;OptionX.Duration=15;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Module Color_IP',OptionX)
% % Add Morphine Significance
% h.NodeCData = g.Nodes.P_adj;
% set(gcf,'Colormap',flipud(inferno),'Color','k');
% caxis([0 .5]);
% % Set up recording parameters (optional), and record
% OptionX.FrameRate=30;OptionX.Duration=20;OptionX.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'MinSpan Morphine Significance_IP',OptionX)

%%
w.loadDESEQ("..\DESeq2\2 DESeq2 EtOH0 IP vs All Air IP.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
%
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
% g = w.removeDisconnectedNodes(g);
mnst = minspantree(g,"Type","forest");

closeness=centrality(g,'closeness','Cost',g.Edges.Weight);
closeness_mnst=centrality(mnst,'closeness','Cost',mnst.Edges.Weight);
degree=centrality(g,'degree');
degree_mnst=centrality(mnst,'degree');

g.Nodes.Closeness=closeness;
g.Nodes.Degree = degree;
mnst.Nodes.Closeness=closeness_mnst;
mnst.Nodes.Degree = degree_mnst;

% % Save Table of Module of Interest
% gTable=sortrows(g.Nodes,'Closeness','descend');
% writetable(gTable,[Interest{i} 'graphSort_FentFive_vs_FentOne.xlsx']);
% for k=1:5
% gHubs{i,k}=gTable.Name{k};
% end

mnstTable=sortrows(mnst.Nodes,'Degree','descend');
writetable(mnstTable,[Interest{i} 'mnstSort_ETOH.xlsx']);
for k=1:5
mnstHubs{i,k}=mnstTable.Name{k};
end

% Centrality
f1=figure('color','w','position',[100 100 300 270]);
scatter(g.Nodes.Closeness,(g.Nodes.Wald_Stats),6,'filled','o');
p = polyfit(g.Nodes.Closeness,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,g.Nodes.Closeness);
hold on
scatter(g.Nodes.Closeness,y1,4,'filled','o');
ylabel('ETOH (Wald Stat)');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=g.Nodes.Closeness(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq,'%.2f') ', P= ' num2str(Pval,'%.2f')]);
savefig(gcf,[Interest{i} ' Centrality ETOH.fig']);
export_fig([Interest{i} ' Centrality ETOH.png'],'-m5');

% 
OptionZ.EdgeAlpha=1;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=18
OptionZ.SigNode=1;
OptionZ.Layout='tree';
h = w.plotGraph(g,OptionZ);
h.NodeCData = mnst.Nodes.Wald_Stats;
h.EdgeColor = [0.9,0.6,0.10];
set(gcf,'Colormap',crameri('Vanimo'),'Color','k','Position',[475,140,500,320]);
c=colorbar('Color','w');
c.Label.String="Wald-Statitic";
c.Label.FontSize=14;
c.FontSize=14;
caxis([-8 8]);
savefig(gcf,[Interest{i} ' Network Tree.fig']);
export_fig([Interest{i} ' Network Tree.png'],'-m5');

% Only Plot 100 Genes
if height(g.Nodes)>30
[~, NodesToRemove] = mink(abs(g.Nodes.Closeness), height(g.Nodes)-30);
g = rmnode(g, NodesToRemove);
end

OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=12;
OptionZ.SigNode=1;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
% h.NodeCData = mnst.Nodes.Wald_Stats;
savefig(gcf,[Interest{i} ' Network Circle.fig']);
export_fig([Interest{i} ' Network Circle.png'],'-m5');

close all
end

%%
w.loadDESEQ("..\DESeq2\3 DESeq2 ETOH8 IP vs ETOH0 IP.xlsx");
dsa = w.deseqTable(:,[2,3,4,5,6,7,9]);
statarray = grpstats(dsa,'moduleColor');
Interest=statarray.moduleColor(abs(statarray.mean_Wald_Stats)>1.5);
Interest=cellstr(Interest);
%
for i=1:length(Interest)
w.createGraphTOM;
g = w.getGraphOfModule(Interest{i,1})
% g = w.removeDisconnectedNodes(g);
mnst = minspantree(g,"Type","forest");

closeness=centrality(g,'closeness','Cost',g.Edges.Weight);
closeness_mnst=centrality(mnst,'closeness','Cost',mnst.Edges.Weight);
degree=centrality(g,'degree');
degree_mnst=centrality(mnst,'degree');

g.Nodes.Closeness=closeness;
g.Nodes.Degree = degree;
mnst.Nodes.Closeness=closeness_mnst;
mnst.Nodes.Degree = degree_mnst;

% % Save Table of Module of Interest
% gTable=sortrows(g.Nodes,'Closeness','descend');
% writetable(gTable,[Interest{i} 'graphSort_FentFive_vs_FentOne.xlsx']);
% for k=1:5
% gHubs{i,k}=gTable.Name{k};
% end

mnstTable=sortrows(mnst.Nodes,'Degree','descend');
writetable(mnstTable,[Interest{i} 'mnstSort_Withdrawal.xlsx']);
for k=1:5
mnstHubs{i,k}=mnstTable.Name{k};
end

% Centrality
f1=figure('color','w','position',[100 100 300 270]);
scatter(g.Nodes.Closeness,(g.Nodes.Wald_Stats),6,'filled','o');
p = polyfit(g.Nodes.Closeness,(g.Nodes.Wald_Stats),1);
y1 = polyval(p,g.Nodes.Closeness);
hold on
scatter(g.Nodes.Closeness,y1,4,'filled','o');
ylabel('Withdrawal (Wald Stat)');
xlabel('Centrality');
%xlim([.001 .003]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
x=g.Nodes.Closeness(~isnan(g.Nodes.Wald_Stats));
y=g.Nodes.Wald_Stats(~isnan(g.Nodes.Wald_Stats));
[Rho Pval]=corr(x,y,'Type','Spearman');
R_sq=Rho^2;
title(['R^2= ' num2str(R_sq,'%.2f') ', P= ' num2str(Pval,'%.2f')]);
savefig(gcf,[Interest{i} ' Centrality Withdrawal.fig']);
export_fig([Interest{i} ' Centrality Withdrawal.png'],'-m5');

% 
OptionZ.EdgeAlpha=1;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=18
OptionZ.SigNode=1;
OptionZ.Layout='tree';
h = w.plotGraph(g,OptionZ);
h.NodeCData = mnst.Nodes.Wald_Stats;
h.EdgeColor = [0.9,0.6,0.10];
set(gcf,'Colormap',crameri('Vanimo'),'Color','k','Position',[475,140,500,320]);
c=colorbar('Color','w');
c.Label.String="Wald-Statitic";
c.Label.FontSize=14;
c.FontSize=14;
caxis([-8 8]);
savefig(gcf,[Interest{i} ' Network Tree.fig']);
export_fig([Interest{i} ' Network Tree.png'],'-m5');

% Only Plot 100 Genes
if height(g.Nodes)>30
[~, NodesToRemove] = mink(abs(g.Nodes.Closeness), height(g.Nodes)-30);
g = rmnode(g, NodesToRemove);
end

OptionZ.EdgeAlpha=.5;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=12;
OptionZ.SigNode=1;
OptionZ.Layout='circle';
h = w.plotGraph(g,OptionZ);
% h.NodeCData = mnst.Nodes.Wald_Stats;
savefig(gcf,[Interest{i} ' Network Circle.fig']);
export_fig([Interest{i} ' Network Circle.png'],'-m5');

close all
end


%% Plot graph of individual gene network
w.createGraphTOM;

for i=1:height(Interest)
for k=1:width(gHubs)    
g = w.getGraphOfNeighbors(gHubs{i,k}, 50);
center=centrality(g,'degree');
g.Nodes.Closeness=center;
% Save Table of Module of Interest
gTable=sortrows(g.Nodes,'Closeness','descend');
writetable(gTable,[Interest{i} '-' gHubs{i,k} '-localGraph.xlsx']);
end
end

for i=1:height(Interest)
for k=1:width(mnstHubs)    
g = w.getGraphOfNeighbors(mnstHubs{i,k}, 50);
center=centrality(g,'degree');
g.Nodes.Closeness=center;
% Save Table of Module of Interest
gTable=sortrows(g.Nodes,'Closeness','descend');
writetable(gTable,[Interest{i} '-' mnstHubs{i,k} '-localMST.xlsx']);
end
end

%% NLRP3
g = w.getGraphOfNeighbors('Nlrp3', 50);
center=centrality(g,'degree');
g.Nodes.Closeness=center;
% Save Table of Module of Interest
gTable=sortrows(g.Nodes,'Closeness','descend');
writetable(gTable,['Nlrp3-localGraph.xlsx']);
OptionZ.EdgeAlpha=1;
OptionZ.LineWidth=1;
OptionZ.NodeFontSize=18
OptionZ.SigNode=1;
OptionZ.Layout='force';
h = w.plotGraph(g,OptionZ);
h.NodeCData = mnst.Nodes.Wald_Stats;
h.EdgeColor = [0.4,0.05,0.10];
set(gcf,'Colormap',crameri('Vanimo'),'Color','k');
c=colorbar('Color','w');
c.Label.String="Wald-Statitic";
c.Label.FontSize=14;
c.FontSize=14;
caxis([-8 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Enrichment Sorting - Doesnt work Yet
d=dir('\\10.159.50.8\kcoffey\Coffey.Lab\Content.Production\Collaborator.Content\2025.04.01 Bergkamp WGCNA Fresh\Raw Figures\*.xlsx')
% Manually Delete Rows LOL...
for i=1:2:(height(d)-1) 
Fent=readtable(fullfile(d(i).folder,d(i).name));
Enrich=readtable(fullfile(d(i+1).folder,d(i+1).name));
% Manually fix Duplicate Probes
zF=zscore(Fent.Wald_Stats);
zE=zscore(Enrich.Wald_Stats);
zZ=zF+zE;
[~, Zi] =sort(zZ,'descend');
probes=Fent.Name(Zi);
p250=probes(1:250);
c=strsplit(d(i).name,'_');
writetable(table(p250),[c{1,1} '_EnrichSort250.xlsx']);
end