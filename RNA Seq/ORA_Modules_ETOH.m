%% 
%% ORA for Orange Yellow
clear all
fig_fold='.\Raw Figures';
AP=readtable(".\Module ORA\Orange Yellow Withdrawal\enriched_geneset_ap_clusters_wg_result1585269850.txt",'Delimiter','\t');
t = readtable(".\Module ORA\Orange Yellow Withdrawal\enrichment_results_wg_result1585269850.txt");
[c2,ia,ib]=intersect(t.geneSet,AP{:,1},'stable');

t=t(ia,:);
t=sortrows(t,'enrichmentRatio','descend');
t.database=categorical(t.database);

if height(t)>16
idx=extractTopIndices(t, 16);
end
t=t(idx,:);
t=sortrows(t,'enrichmentRatio','descend');
t.geneSet=regexprep(t.geneSet, '_', ' ');
t.description(cellfun(@isempty,t.description))=t.geneSet(cellfun(@isempty,t.description));

for i=1:height(t)
    cellContents = t.description{i};
    if length(cellContents)~=0
        if length(cellContents)>35
            t.description{i} = [cellContents(1:35)  ' (' num2str(t.overlap(i)) ')'];
        else
            t.description{i} = [cellContents(1:end) ' (' num2str(t.overlap(i)) ')'];
        end
    end
end

t.geneSet=categorical(t.geneSet);
t.description=categorical(t.description);

% Plotting Defferntial Expression
[~, rgb] = colornames('Crayola','Orange Yellow');
f1=figure('color','w','position',[100 100 350 350]);
b=barh(flipud(t.enrichmentRatio));
yticks([1:height(t)]);
yticklabels(flipud(t.description));
box on;
set(gca,'TickDir','none','LineWidth',1.5,'YAxisLocation','right');
b.BarWidth=.9;
b.FaceColor=rgb;
b.EdgeColor=[0,0,0];
xlabel('Enrichment Ratio');
savefig(gcf,fullfile(fig_fold, 'Orange Yellow ORA.fig'));
export_fig(fullfile(fig_fold, 'Orange Yellow ORA.png'), '-m3'); % Save the Figure

%% ORA for Outer Space
clear all
clear all
fig_fold='.\Raw Figures';
AP=readtable(".\Module ORA\Outer Space Withdrawal\enriched_geneset_ap_clusters_wg_result1585269817.txt",'Delimiter','\t');
t = readtable(".\Module ORA\Outer Space Withdrawal\enrichment_results_wg_result1585269817.txt");
[c2,ia,ib]=intersect(t.geneSet,AP{:,1},'stable');

t=t(ia,:);
t=sortrows(t,'enrichmentRatio','descend');
t.database=categorical(t.database);

if height(t)>16
    idx=extractTopIndices(t, 16);
end
t=t(idx,:);
t=sortrows(t,'enrichmentRatio','descend');
t.geneSet=regexprep(t.geneSet, '_', ' ');
t.description(cellfun(@isempty,t.description))=t.geneSet(cellfun(@isempty,t.description));

for i=1:height(t)
    cellContents = t.description{i};
    if length(cellContents)~=0
        if length(cellContents)>35
            t.description{i} = [cellContents(1:35)  ' (' num2str(t.overlap(i)) ')'];
        else
            t.description{i} = [cellContents(1:end) ' (' num2str(t.overlap(i)) ')'];
        end
    end
end

t.geneSet=categorical(t.geneSet);
t.description=categorical(t.description);

% Plotting Defferntial Expression
[~, rgb] = colornames('Crayola','Outer Space');
f1=figure('color','w','position',[100 100 350 350]);
b=barh(flipud(t.enrichmentRatio));
yticks([1:height(t)]);
yticklabels(flipud(t.description));
box on;
set(gca,'TickDir','none','LineWidth',1.5,'YAxisLocation','right');
b.BarWidth=.9;
b.FaceColor=rgb;
b.EdgeColor=[0,0,0];
xlabel('Enrichment Ratio');
savefig(gcf,fullfile(fig_fold, 'Outer Space ORA.fig'));
export_fig(fullfile(fig_fold, 'Outer Space ORA.png'), '-m3'); % Save the Figure

%% ORA for Dandelion
clear all
clear all
fig_fold='.\Raw Figures';
AP=readtable(".\Module ORA\DandelionORA_2025\enriched_geneset_ap_clusters_wg_result1745179793.txt",'Delimiter','\t');
t = readtable(".\Module ORA\DandelionORA_2025\enrichment_results_wg_result1745179793.txt");
[c2,ia,ib]=intersect(t.geneSet,AP{:,1},'stable');

t=t(ia,:);
t=sortrows(t,'enrichmentRatio','descend');
t.database=categorical(t.database);

if height(t)>16
    idx=extractTopIndices(t, 16);
    t=t(idx,:);
end
t=sortrows(t,'enrichmentRatio','descend');
t.geneSet=regexprep(t.geneSet, '_', ' ');
t.description(cellfun(@isempty,t.description))=t.geneSet(cellfun(@isempty,t.description));

%
for i=1:height(t)
    cellContents = t.description{i};
    if length(cellContents)~=0
        if length(cellContents)>35
            t.description{i} = [cellContents(1:35)  ' (' num2str(t.overlap(i)) ')'];
        else
            t.description{i} = [cellContents(1:end) ' (' num2str(t.overlap(i)) ')'];
        end
    end
end

t.geneSet=categorical(t.geneSet);
t.description=categorical(t.description);

% Plotting Defferntial Expression
[~, rgb] = colornames('Crayola','Dandelion');
f1=figure('color','w','position',[100 100 350 350]);
b=barh(flipud(t.enrichmentRatio));
yticks([1:height(t)]);
yticklabels(flipud(t.description));
box on;
set(gca,'TickDir','none','LineWidth',1.5,'YAxisLocation','right');
b.BarWidth=.9;
b.FaceColor=rgb;
b.EdgeColor=[0,0,0];
xlabel('Enrichment Ratio');
savefig(gcf,fullfile(fig_fold, 'Dandelion ORA.fig'));
export_fig(fullfile(fig_fold, 'Dandelion ORA.png'), '-m3'); % Save the Figure

%% ORA for Turquoise Blue
clear all
clear all
fig_fold='.\Raw Figures';
AP=readtable(".\Module ORA\Turquoise Blue ETOH\enriched_geneset_ap_clusters_wg_result1585269781.txt",'Delimiter','\t');
t = readtable(".\Module ORA\Turquoise Blue ETOH\enrichment_results_wg_result1585269781.txt");
[c2,ia,ib]=intersect(t.geneSet,AP{:,1},'stable');

t=t(ia,:);
t=sortrows(t,'enrichmentRatio','descend');
t.database=categorical(t.database);

if height(t)>16
    idx=extractTopIndices(t, 16);
    t=t(idx,:);
end
t=sortrows(t,'enrichmentRatio','descend');
t.geneSet=regexprep(t.geneSet, '_', ' ');
t.description(cellfun(@isempty,t.description))=t.geneSet(cellfun(@isempty,t.description));

%
for i=1:height(t)
    cellContents = t.description{i};
    if length(cellContents)~=0
        if length(cellContents)>35
            t.description{i} = [cellContents(1:35)  ' (' num2str(t.overlap(i)) ')'];
        else
            t.description{i} = [cellContents(1:end) ' (' num2str(t.overlap(i)) ')'];
        end
    end
end

t.geneSet=categorical(t.geneSet);
t.description=categorical(t.description);

% Plotting Defferntial Expression
[~, rgb] = colornames('Crayola','Turquoise Blue');
f1=figure('color','w','position',[100 100 350 350]);
b=barh(flipud(t.enrichmentRatio));
yticks([1:height(t)]);
yticklabels(flipud(t.description));
box on;
set(gca,'TickDir','none','LineWidth',1.5,'YAxisLocation','right');
b.BarWidth=.9;
b.FaceColor=rgb;
b.EdgeColor=[0,0,0];
xlabel('Enrichment Ratio');
savefig(gcf,fullfile(fig_fold, 'Turquoise Blue ORA.fig'));
export_fig(fullfile(fig_fold, 'Turquoise Blue ORA.png'), '-m3'); % Save the Figure

%%
function idxResult = extractTopIndices(tbl, n)
    % Check input
    if ~istable(tbl) || ~all(ismember({'enrichmentRatio', 'database'}, tbl.Properties.VariableNames))
        error('Input must be a table with variables "enrichmentRatio" and "database".');
    end

    % Get list of unique databases (categories)
    databases = unique(tbl.database);
    numDatabases = numel(databases);
    
    % For each database, sort rows by descending enrichmentRatio and store indices
    sortedIndices = cell(numDatabases, 1);
    for i = 1:numDatabases
        db = databases(i);
        rows = find(tbl.database == db);
        [~, sortOrder] = sort(tbl.enrichmentRatio(rows), 'descend');
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
