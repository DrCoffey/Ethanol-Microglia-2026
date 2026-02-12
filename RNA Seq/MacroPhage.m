t=readtable('Macrophage\Macrophage.xlsx');
t.Gene=categorical(t.Gene);
t.CellType=categorical(t.CellType);

% Import Raw Gene Counts
rt=readtable('DataSheets\Gene Table IP.csv');
rt.Genes=categorical(rt.Genes);

% Raw Counts
for i=1:height(t)
idx=find(rt.Genes==t.Gene(i));
mCount(i,1)=mean(mean(rt{idx,2:end}));
end
t=[t table(mCount)];

% Enrichment
nt=readtable('DeSeq2\4 DESeq2 All IP vs All IN.xlsx');
nt.GeneID=categorical(nt.GeneID);
for i=1:height(t)
idx=find(nt.GeneID==t.Gene(i));
Wald(i,1)=nt.Wald_Stats(idx);
pVal(i,1)=nt.P_adj(idx);
end
t=[t table(Wald, pVal)];

f1=figure(Position=[1200 500 850 400]);
g=gramm('x',t.Gene,'y',(t.mCount+1),'color',t.CellType);
g.geom_bar(width=.9);
g.set_names('x','Gene','y','Mean TPM','color','Cell Type');
g.set_order_options('x',t.Gene);
g.axe_property('YLim',[0 10000],'TickDir','out','LineWidth',1.5,'FontSize',12)
g.draw();
exportgraphics(f1,'Macrophage/macroCount.png','Resolution',600);

f2=figure(Position=[1200 500 400 400]);
g=gramm('x',t.CellType,'y',t.Wald,'color',t.CellType);
g.stat_boxplot();
g.set_names('x','Cell Type','y','Enrichment (Wald)','color','Cell Type');
g.axe_property('YLim',[-2 10],'TickDir','out','LineWidth',1.5,'FontSize',12,'XLabel',[])
g.set_order_options('x',{'Microglia', 'Macrophage'});
g.draw();
exportgraphics(f2,'Macrophage/macroWald.png','Resolution',600);

f3=figure(Position=[1200 500 400 400]);
g=gramm('x',t.CellType,'y',t.pVal,'color',t.CellType);
g.stat_boxplot();
g.set_names('x','Cell Type','y','Enrichment log(qValue)','color','Cell Type');
g.axe_property('YLim',[0 1],'TickDir','out','LineWidth',1.5,'FontSize',12,'YScale','log','XLabel',[])
g.set_order_options('x',{'Microglia', 'Macrophage'});
g.draw();
hold on
yline(g.facet_axes_handles,.1,'--','q < .1','LineWidth',1.5);
exportgraphics(f3,'Macrophage/macroP.png','Resolution',600);