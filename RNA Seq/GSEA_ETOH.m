%% 
%% GSEA Enrichment
clear all
fig_fold='.\Raw Figures';
t = readtable(".\GSEA\Enrichment\Molecular Function Enrichment.xlsx");
t.NES_ABS=abs(t.NES);

t=sortrows(t,'NES_ABS','descend');

for i=1:height(t)
    cellContents = t.Description{i};
    if length(cellContents)~=0
        if length(cellContents)>35
            t.Description{i} = [cellContents(1:35)  ' (' num2str(t.LeadingEdgeNumber(i)) ')'];
        else
            t.Description{i} = [cellContents(1:end) ' (' num2str(t.LeadingEdgeNumber(i)) ')'];
        end
    end
end

t.GeneSet=categorical(t.GeneSet);
t.Description=categorical(t.Description);

t=t(1:24,:);
t=sortrows(t,'NES','descend');

for i=1:height(t)
    if t.NES(i)>0
       cvars(i)=.05-t.FDR(i);
    else
       cvars(i)=-.05+t.FDR(i); 
    end
end

% Plotting Defferntial Expression
[~, rgb] = colornames('Crayola','Orange');
f1=figure('color','w','position',[100 100 500 400]);
b=barh(flipud(t.NES));
yticks([1:height(t)]);
yticklabels(flipud(t.Description));
box on;
set(gca,'TickDir','none','LineWidth',1,'YAxisLocation','right');
b.BarWidth=.9;
cmap=crameri('Vanimo');
cvars_map = (cvars - min(cvars))/(max(cvars)-min(cvars)) * (256-1) + 1;
cdat=cmap(floor(cvars_map), :);
b.FaceColor = 'flat';
b.CData = flipud(cdat);
%b.FaceColor=rgb;
b.EdgeColor=[1,1,1];
xlabel('Enrichment Ratio');
%colorbar;
savefig(gcf,fullfile(fig_fold, 'Enrichment GSEA.fig'));
export_fig(fullfile(fig_fold, 'Enrichment GSEA.png'), '-m3'); % Save the Figure

hf = figure('Units','normalized'); 
colormap(cmap)
hCB = colorbar('north');
set(gca,'Visible',false)
hCB.Position = [0.15 0.3 0.74 0.4];
hf.Position(4) = 0.1000;
export_fig(fullfile(fig_fold, 'ColorBar.png'), '-m3'); % Save the Figure

%% ETOH GSEA with Affinity Propagation
clear all
fig_fold='.\Raw Figures';
AP=readtable(".\GSEA\IP\Ethanol v Air\Mol Fun\enriched_geneset_ap_clusters_wg_result1568400084.txt",'Delimiter','\t');
t = readtable(".\GSEA\IP\Ethanol v Air\Mol Fun\enrichment_results_wg_result1568400084.txt");
[c2,ia,ib]=intersect(t.geneSet,AP{:,1},'stable');
t=t(ia,:);
database=repmat({'GO_MF'},[height(t),1]);
t.database=categorical(database);

AP=readtable(".\GSEA\IP\Ethanol v Air\Reactome\enriched_geneset_ap_clusters_wg_result1568400113.txt",'Delimiter','\t');
t2 = readtable(".\GSEA\IP\Ethanol v Air\Reactome\enrichment_results_wg_result1568400113.txt");
[c2,ia,ib]=intersect(t2.geneSet,AP{:,1},'stable');
t2=t2(ia,:);
database=repmat({'Reactome'},[height(t2),1]);
t2.database=categorical(database);

AP=readtable(".\GSEA\IP\Ethanol v Air\Wiki\enriched_geneset_ap_clusters_wg_result1586971500.txt",'Delimiter','\t');
t3 = readtable(".\GSEA\IP\Ethanol v Air\Wiki\enrichment_results_wg_result1586971500.txt");
[c2,ia,ib]=intersect(t3.geneSet,AP{:,1},'stable');
t3=t3(ia,:);
database=repmat({'Wiki'},[height(t3),1]);
t3.database=categorical(database);

t=[t;t2;t3];

t.NES_ABS=abs(t.normalizedEnrichmentScore);
t=sortrows(t,'NES_ABS','descend');
te=t(t.normalizedEnrichmentScore>0,:);
td=t(t.normalizedEnrichmentScore<0,:);

if height(t)>16
idxE=extractTopIndices(te, 12);
idxD=extractTopIndices(td, 12);
end

t=[te(idxE,:);td(idxD,:)];
t=sortrows(t,'normalizedEnrichmentScore','descend');

for i=1:height(t)
    cellContents = t.description{i};
    if length(cellContents)~=0
        if length(cellContents)>45
            t.description{i} = [cellContents(1:45)  ' (' num2str(t.leadingEdgeNum(i)) ')'];
        else
            t.description{i} = [cellContents(1:end) ' (' num2str(t.leadingEdgeNum(i)) ')'];
        end
    end
end

t.geneSet=categorical(t.geneSet);
t.description=categorical(t.description);

for i=1:height(t)
    if t.normalizedEnrichmentScore(i)>0
       cvars(i,1)=max(t.FDR)-t.FDR(i);
    else
       cvars(i,1)=-max(t.FDR)+t.FDR(i); 
    end
end

% Plotting Defferntial Expression
f1=figure('color','w','position',[100 100 650 400]);
b=barh(flipud(t.normalizedEnrichmentScore));
yticks([1:height(t)]);
yticklabels(flipud(t.description));
box on;
set(gca,'TickDir','none','LineWidth',1,'YAxisLocation','right');
b.BarWidth=.9;
cmap=crameri('Vanimo');
cvars_map = (cvars - min(cvars))/(max(cvars)-min(cvars)) * (256-1) + 1;
cdat=cmap(floor(cvars_map), :);
b.FaceColor = 'flat';
b.CData = flipud(cdat);
b.EdgeColor=[1,1,1];
xlabel('Enrichment Ratio');
%colorbar;
savefig(gcf,fullfile(fig_fold, 'ETOH GSEA.fig'));
export_fig(fullfile(fig_fold, 'ETOH GSEA.png'), '-m3'); % Save the Figure


%% Withdrawal GSEA with Affinity Propagation
clear all
fig_fold='.\Raw Figures';
AP=readtable(".\GSEA\IP\Withdrawal v Ethanol\Mol Fun\enriched_geneset_ap_clusters_wg_result1568401775.txt",'Delimiter','\t');
t = readtable(".\GSEA\IP\Withdrawal v Ethanol\Mol Fun\enrichment_results_wg_result1568401775.txt");
[c2,ia,ib]=intersect(t.geneSet,AP{:,1},'stable');
t=t(ia,:);
database=repmat({'GO_MF'},[height(t),1]);
t.database=categorical(database);

AP=readtable(".\GSEA\IP\Withdrawal v Ethanol\Reactome\enriched_geneset_ap_clusters_wg_result1568401789.txt",'Delimiter','\t');
t2 = readtable(".\GSEA\IP\Withdrawal v Ethanol\Reactome\enrichment_results_wg_result1568401789.txt");
[c2,ia,ib]=intersect(t2.geneSet,AP{:,1},'stable');
t2=t2(ia,:);
database=repmat({'Reactome'},[height(t2),1]);
t2.database=categorical(database);

AP=readtable(".\GSEA\IP\Withdrawal v Ethanol\Wiki\enriched_geneset_ap_clusters_wg_result1568401795.txt",'Delimiter','\t');
t3 = readtable(".\GSEA\IP\Withdrawal v Ethanol\Wiki\enrichment_results_wg_result1568401795.txt",'Delimiter','\t');
[c2,ia,ib]=intersect(t3.geneSet,AP{:,1},'stable');
t3=t3(ia,:);
database=repmat({'Wiki'},[height(t3),1]);
t3.database=categorical(database);

t=[t;t2;t3];

t.NES_ABS=abs(t.normalizedEnrichmentScore);
t=sortrows(t,'NES_ABS','descend');
te=t(t.normalizedEnrichmentScore>0,:);
td=t(t.normalizedEnrichmentScore<0,:);

if height(t)>16
idxE=extractTopIndices(te, 12);
idxD=extractTopIndices(td, 12);
end

t=[te(idxE,:);td(idxD,:)];
t=sortrows(t,'normalizedEnrichmentScore','descend');

for i=1:height(t)
    cellContents = t.description{i};
    if length(cellContents)~=0
        if length(cellContents)>45
            t.description{i} = [cellContents(1:45)  ' (' num2str(t.leadingEdgeNum(i)) ')'];
        else
            t.description{i} = [cellContents(1:end) ' (' num2str(t.leadingEdgeNum(i)) ')'];
        end
    end
end
t.geneSet=categorical(t.geneSet);
t.description=categorical(t.description);

for i=1:height(t)
    if t.normalizedEnrichmentScore(i)>0
       cvars(i,1)=max(t.FDR)-t.FDR(i);
    else
       cvars(i,1)=-max(t.FDR)+t.FDR(i); 
    end
end

% Plotting Defferntial Expression
f1=figure('color','w','position',[100 100 650 400]);
b=barh(flipud(t.normalizedEnrichmentScore));
yticks([1:height(t)]);
yticklabels(flipud(t.description));
box on;
set(gca,'TickDir','none','LineWidth',1,'YAxisLocation','right');
b.BarWidth=.9;
cmap=crameri('Vanimo');
cvars_map = (cvars - min(cvars))/(max(cvars)-min(cvars)) * (256-1) + 1;
cdat=cmap(floor(cvars_map), :);
b.FaceColor = 'flat';
b.CData = flipud(cdat);
b.EdgeColor=[1,1,1];
xlabel('Enrichment Ratio');
%colorbar;
savefig(gcf,fullfile(fig_fold, 'Withdrawal GSEA.fig'));
export_fig(fullfile(fig_fold, 'Withdrawal GSEA.png'), '-m3'); % Save the Figure


%%
function idxResult = extractTopIndices(tbl, n)
    % Check input
    if ~istable(tbl) || ~all(ismember({'NES_ABS', 'database'}, tbl.Properties.VariableNames))
        error('Input must be a table with variables "NES_ABS" and "database".');
    end

    % Get list of unique databases (categories)
    databases = unique(tbl.database);
    numDatabases = numel(databases);
    
    % For each database, sort rows by descending NES_ABS and store indices
    sortedIndices = cell(numDatabases, 1);
    for i = 1:numDatabases
        db = databases(i);
        rows = find(tbl.database == db);
        [~, sortOrder] = sort(tbl.NES_ABS(rows), 'descend');
        sortedIndices{i} = rows(sortOrder); % Store sorted indices
    end

    % Initialize tracking of current position for each database (start at 1)
    currentPos = ones(numDatabases, 1); % Tracks next index for each database (starting at 1)
    idxResult = [];
    totalExtracted = 0;

    % Round-robin extraction
    while totalExtracted < n
        allExhausted = true;
        for i = 1:numDatabases
            if totalExtracted >= n
                break;
            end
            if currentPos(i) <= numel(sortedIndices{i})
                idxResult(end+1, 1) = sortedIndices{i}(currentPos(i)); %#ok<AGROW>
                currentPos(i) = currentPos(i) + 1;
                totalExtracted = totalExtracted + 1;
                allExhausted = false;
            end
        end
        if allExhausted
            break; % no more values to extract
        end
    end
end
