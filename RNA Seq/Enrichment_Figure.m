% Sets Directory For Data & Figures
clear all 
close all
fig_fold='.\Raw Figures';
t=readtable(".\DESeq2\4 DESeq2 All IP vs All IN.xlsx"); % Import seq data

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 350 310]);
hold on;
scatter(log2(t.BaseMean),t.log2_FC_,6,[.5 .5 .5],'filled','o');
cvars=t.Wald_Stats;
cvars_map = (cvars - min(cvars))/(max(cvars)-min(cvars)) * (256-1) + 1;
cmap=crameri('Vanimo');
idx=t.P_adj<.10;
cdat=cmap(floor(cvars_map(idx)), :);
scatter(log2(t.BaseMean(idx)),t.log2_FC_(idx),8,cdat,'filled','o','markerfacealpha',.6);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
ylim([-6 8]);
xlim([-5 15]);
colormap(cmap);
h = colorbar;
set(get(h,'title'),'string','(Wald)');
% Create ylabel
ylabel('log2(Fold Change)');
xlabel('log2(TPM)');

h.TickLabels=floor(y);
t2=sortrows(t,'Wald_Stats','descend');
t2=t2(~isnan(t2.Wald_Stats),:);
t2=t2(1:20,:);
s2=scatter(log2(t2.BaseMean),t2.log2_FC_,(35:-1:16),cdat(1:20,:),'filled','o','markerfacealpha',1,'MarkerEdgeColor','k');
sgene=t2.GeneID(1:20);
savefig(gcf,fullfile(fig_fold, 'Enrichment DESeq.fig'));
export_fig(fullfile(fig_fold, 'Enrichment DESeq.png'), '-m3'); % Save the Figure

f1=figure('color','w','position',[100 100 350 310]);
s2=scatter(ones([20,1]),(20:-1:1)',(35:-1:16),cdat(1:20,:),'filled','o','markerfacealpha',1,'MarkerEdgeColor','k');
text(ones([20,1])+.2,(20:-1:1)',sgene);
savefig(gcf,fullfile(fig_fold, 'Enrichment DESeq leg.fig'));
export_fig(fullfile(fig_fold, 'Enrichment DESeq leg.png'), '-m3'); % Save the Figure

%% Between Experiment Comparison
opts=detectImportOptions(".\DESeq2\4 DESeq2 All IP vs All IN.xlsx")
etoh=readtable(".\DESeq2\4 DESeq2 All IP vs All IN.xlsx",opts);
morph=readtable(".\DESeq2\DESeq2_IP_Input - Morphine.xlsx",opts);

[c ia ib]=intersect(morph.GeneID,etoh.GeneID);
% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 300 270]);
hold on;
plot([-6 8],[-6 8],'--k');
scatter((morph.log2_FC_(ia)),(etoh.log2_FC_(ib)),6,'o','filled');
box off;
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');

% Individual Significance
msig=morph(morph.P_adj<.05 & morph.log2_FC_>0,:);
esig=etoh(etoh.P_adj<.05 & etoh.log2_FC_>0,:);

[c ia ib]=intersect(msig.GeneID,etoh.GeneID);
scatter((msig.log2_FC_(ia)),(etoh.log2_FC_(ib)),6,'o','filled');

[c ia ib]=intersect(morph.GeneID,esig.GeneID);
scatter((morph.log2_FC_(ia)),(esig.log2_FC_(ib)),6,'o','filled');

[c ia ib]=intersect(msig.GeneID,esig.GeneID);
scatter((msig.log2_FC_(ia)),(esig.log2_FC_(ib)),6,'o','filled');

% legend({'Ref Line','All Genes','Morphine Enrichment','ETOH Enrichment','Both'},'Location','eastoutside');
xlabel('Morphine Enrichment log2(FC)', 'FontSize',12);
ylabel('ETOH Enrichment log2(FC)','FontSize',12);

savefig(gcf,fullfile(fig_fold, 'Enrichment Across Experiments.fig'));
export_fig(fullfile(fig_fold, 'Enrichment Across Experiments.png'), '-m3'); % Save the Figure

both_sig=esig(ib,:);
c=0;
for i=0:.5:6.5
    c=c+1;
    propsig(c)=sum(both_sig.log2_FC_ > i & both_sig.log2_FC_ < i+.5)/sum(esig.log2_FC_ > i & esig.log2_FC_ < i+.5)
end

f1=figure('color','w','position',[100 100 300 270]);
box off;
hold on;
p=plot([.5:.5:7],propsig*100,'-g','LineWidth',2);
plot([.5,7],[(height(both_sig)/height(esig)*100) (height(both_sig)/height(esig)*100)],'--k','LineWidth',1);
ylim([0 100]);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
p.Color=[0.4660 0.6740 0.1880];
xlabel('ETOH Enrichment log2(FC)', 'FontSize',12);
ylabel('Morphine Concordance (%)','FontSize',12);

savefig(gcf,fullfile(fig_fold, 'Enrichment Concordance.fig'));
export_fig(fullfile(fig_fold, 'Enrichment Concordance.png'), '-m3'); % Save the Figure
