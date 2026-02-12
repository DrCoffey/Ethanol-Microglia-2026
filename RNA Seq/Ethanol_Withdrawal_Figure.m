% Sets Directory For Data & Figures
clear all 
close all
fig_fold='.\Raw Figures';
files=dir('.\DESeq2\*DESeq2*'); %find FST data files

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 1100 280]);
ha = tight_subplot(1,3,[.075 .075],[.2 .1],[.075 .075]) ;
% maxW=max(t.Wald_Stats);
% minW=min(t.Wald_Stats);
maxW=10;
minW=-10; % Better graphing, couple outliers at -17

for l=1:3
t=readtable(fullfile(files(l).folder,files(l).name)); % Import seq data    
axes(ha(l));
scatter(log2(t.BaseMean),t.log2_FC_,6,[.5 .5 .5],'filled','o');
cvars=t.Wald_Stats;
cvars_map = (cvars - minW)/(maxW-minW) * (256-1) + 1;
cvars_map(cvars_map<1)=1;
cvars_map(cvars_map>256)=256;
cmap=crameri('Vanimo');
idx=t.P_adj<.10;
cdat=cmap(floor(cvars_map(idx)), :);
hold on
scatter(log2(t.BaseMean(idx)),t.log2_FC_(idx),8,cdat,'filled','o','markerfacealpha',.6);
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
ylim([-4 4]);
xlim([-5 15]);
colormap(cmap);
h = colorbar;
set(get(h,'title'),'string','(Wald)');
y=linspace(minW,maxW,6)
h.TickLabels=floor(y);

if l==1
ylabel('log_2(Fold Change)');
text(-4,4,'Air 0h vs Air 8h','FontSize',12)
% Create xlabel
xlabel('log2(TPM)');
end


if l==2
text(-4,4,'Ethanol 0h vs All Air','FontSize',12)
xlabel('log2(TPM)');

end

if l==3
text(-4,4,'Ethanol 8h vs Ethanol 0h','FontSize',12)
xlabel('log2(TPM)');
end
end

savefig(gcf,fullfile(fig_fold, 'ETOH Withdrawal DESeq.fig'));
export_fig(fullfile(fig_fold, 'ETOH Withdrawal DESeq.png'), '-m3'); % Save the Figure


%% IP Morphine v Withdrawal
tmp1=importdata(fullfile(files(2).folder,files(2).name)); % Import seq data
filename1=strtok(files(2).name,'.')    
C1 = strsplit(filename1,'_');
tmp1.textdata=tmp1.textdata(2:end,1);
[a i]=sort(tmp1.textdata);

tmp2=importdata(fullfile(files(3).folder,files(3).name)); % Import seq data
filename2=strtok(files(3).name,'.')    
C2 = strsplit(filename2,'_');
tmp2.textdata=tmp2.textdata(2:end,1);
[b i2]=sort(tmp2.textdata);

% Plotting Defferntial Expression 
f1=figure('color','w','position',[100 100 280 270]);
scatter(tmp1.data(i,2),tmp2.data(i2,2),8,[.5 .5 .5],'filled','o');
set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');

x=tmp1.data(i,2);
y=tmp2.data(i2,2);
c=mean([tmp1.data(i,6) tmp2.data(i2,6)]')';
idx=(c<.10 );
hold on;
scatter(x(idx),y(idx),8,c(idx),'filled','o');
ylabel('Withdrawal log2(FC)');
xlabel('Ethanol log2(FC)');
xlim([-4 4])
ylim([-4 4])
%l=lsline;
colormap(flipud(cmap))
caxis([0 .1])

plot([-4 4],[0 0],'--k');
plot([0 0],[-4 4],'--k');
% r_squared = corr(x(idx),y(idx))^2
% text(-3.8,3.8,strcat("Sig. Gene R^2 = ",num2str(r_squared)))
% r_squared = corr(tmp1.data(i,2),tmp2.data(i2,2))^2
% % text(-3.8,3.4,strcat("All Gene R^2 = ",num2str(r_squared)))
% l(2).Visible="off";

savefig(gcf,fullfile(fig_fold, 'ETOH vs Withdrawal Corr.fig'));
export_fig(fullfile(fig_fold, 'ETOH vs Withdrawal Corr.png'), '-m3'); % Save the Figure


%% for Input
% tmp1=importdata(files(5).name); % Import seq data
% filename1=strtok(files(5).name,'.')    
% C1 = strsplit(filename,'_');
% tmp1.textdata=tmp1.textdata(2:end,1);
% [a i]=sort(tmp1.textdata);
% 
% tmp2=importdata(files(6).name); % Import seq data
% filename2=strtok(files(6).name,'.')    
% C2 = strsplit(filename,'_');
% tmp2.textdata=tmp2.textdata(2:end,1);
% [b i2]=sort(tmp2.textdata);
% 
% % Plotting Defferntial Expression 
% f1=figure('color','w','position',[100 100 450 400]);
% scatter(tmp1.data(i,2),tmp2.data(i2,2),8,[.5 .5 .5],'filled','o');
% set(gca,'FontSize',12,'LineWidth',1.5,'TickDir','out');
% 
% x=tmp1.data(i,2);
% y=tmp2.data(i2,2);
% c=mean([tmp1.data(i,6) tmp2.data(i2,6)]')';
% idx=(c<.10 );
% hold on;
% scatter(x(idx),y(idx),8,c(idx),'filled','o');
% ylabel('Withdrawal IN log_2(Fold Change)');
% xlabel('Morphine Escalation IN log_2(Fold Change)');
% xlim([-4 4])
% ylim([-4 4])
% l=lsline;
% colormap(flipud(cmap))
% caxis([0 .1])
% 
% plot([-4 4],[0 0],'--k');
% plot([0 0],[-4 4],'--k');
% r_squared = corr(x(idx),y(idx))^2
% text(-3.8,3.8,strcat("Sig. Gene R^2 = ",num2str(r_squared)))
% l(2).Visible="off";
% 
% pngFileName = 'Morphine v Withdrawal Input.png'; % Set the File name 
% fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
% export_fig(fullFileName, '-m5'); % Save the Figure
% 
% 
% 
% % Plotting colorbar Expression 
% f2=figure('color','w','position',[100 100 300 800]);
% set(gca,'FontSize',16,'FontWeight','bold','LineWidth',2,'TickDir','out');
% h = colorbar;
% colormap(flipud(jet))
% set(get(h,'title'),'string','(FDR)');
% pngFileName = 'Colorbar.png'; % Set the File name 
% fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
% export_fig(fullFileName, '-m5'); % Save the Figure
% 
% % Piechart
% f1=figure('color','w','position',[100 100 1000 600]);
% set(gca,'FontSize',16,'FontWeight','bold');
% for l=1:4
% tmp=importdata(files(l).name); % Import seq data
% filename=strtok(files(l).name,'.');    
% C = strsplit(filename,'_');
% tmp.textdata=tmp.textdata(2:end,1);  
% if l==1
%    SFG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1);
% elseif l==2
%    SMG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1); 
% elseif l==3
%    UFG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1);  
% elseif l==4
%    UMG=tmp.textdata(tmp.data(:,6)<=.2 & tmp.data(:,2)>0,1);   
% end
% end
% 
% AllGenes=[SFG;SMG;UFG;UFG];
% UniqGenes=unique(AllGenes);
% SFG_i=intersect(unique([SMG;UFG;UFG]),SFG);
% SFG_un=length(SFG)-length(SFG_i);
% SMG_i=intersect(unique([SFG;UFG;UFG]),SMG);
% SMG_un=length(SMG)-length(SMG_i);
% UFG_i=intersect(unique([SMG;SFG;UMG]),UFG);
% UFG_un=length(UFG)-length(UFG_i);
% UMG_i=intersect(unique([SMG;SFG;UFG]),UMG);
% UMG_un=length(UMG)-length(UMG_i);
% sumUnGenes=sum([SFG_un SMG_un UFG_un UMG_un]);
% X=[(length(UniqGenes)-sumUnGenes), SFG_un, SMG_un, UFG_un, UMG_un];
% h = pie(X);
% colormap(jet(5));
% hText = findobj(h,'Type','text'); % text object handles
% percentValues = get(hText,'String'); % percent values
% txt = {'Overlap (1761): ';'Stressed F (252): ';'Stressed M (245): ';'Unstressed F (46): ';'Unstressed M (307): '}; % strings
% combinedtxt = strcat(txt,percentValues); % strings and percent values
% oldExtents_cell = get(hText,'Extent');
% oldExtents = cell2mat(oldExtents_cell); % numeric array
% hText(1).String = combinedtxt(1);
% hText(2).String = combinedtxt(2);
% hText(3).String = combinedtxt(3);
% hText(4).String = combinedtxt(4);
% hText(5).String = combinedtxt(5);
% newExtents_cell = get(hText,'Extent'); % cell array
% newExtents = cell2mat(newExtents_cell); % numeric array 
% width_change = newExtents(:,3)-oldExtents(:,3);
% signValues = sign(oldExtents(:,1));
% offset = signValues.*(width_change/2);
% textPositions_cell = get(hText,{'Position'}); % cell array
% textPositions = cell2mat(textPositions_cell); % numeric array
% textPositions(:,1) = textPositions(:,1) + offset; % add offset 
% hText(1).Position = textPositions(1,:);
% hText(2).Position = textPositions(2,:);
% hText(3).Position = textPositions(3,:);
% hText(4).Position = textPositions(4,:);
% hText(5).Position = textPositions(5,:);
% set(hText,'FontWeight','bold','FontSize',13);
% set(h,'LineWidth',2);
% pngFileName = 'ControlPie.png'; % Set the File name 
% fullFileName = fullfile(fig_fold, pngFileName); % Add Figure Path
% export_fig(fullFileName, '-m5'); % Save the Figure
% 
% xlswrite('SignalGenes.xlsx',UniqGenes);