% common_clustering.m

% requirements: [common_data_folder]/common_data/samplex_data.mat
%                and conneurons.csv
% where x stands for sample number.

% reads samplex_data.mat and makes outperm
% and saves it in [common_cell_order_folder]/outperm.mat 
% several data about A and B neuron groups are also saved in
% [common_cell_order_folder]




%%%%%  Preset variables  %%%%%

samples = 1:24;
project_folder = '/home/iino/gKDR_GMM_50mM/Linked_gKDR_freerun/publish_version';
common_data_folder = fullfile(project_folder, 'common_data'); %'../ver92_Hbase/common_data';
utilities_folder = fullfile(project_folder, 'utilities'); %'../'
metadata_folder = fullfile(project_folder, 'metadata'); %'../'
common_cell_order_folder = fullfile(project_folder, 'common_cell_order'); %'../ver92_Hbase/common_cell_order';

time_step = 5;
method = 'single';


%%%%%  Initial treatment  %%%%%

conNames = table2cell(readtable(fullfile(metadata_folder, 'conneurons.csv'),'ReadVariableNames',false));

if ~exist(common_cell_order_folder, 'dir')
    mkdir(common_cell_order_folder)
end

nsamples = length(samples);
addpath(utilities_folder);

ABclusters = cell(2*nsamples,1);    %% member list for each cluster, two for each sample
ABclustermat = logical(zeros(2*nsamples,length(conNames)));  %% boolean list of members.


%%%%% load data for each sample and clustering %%%%%

disp('<< Simple clustering >>')
for sampleID = samples
    
    disp(['< sample ' num2str(sampleID) ' >'])
    
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'data', 'uniqNames', 'targetcellnames', 'targetcells', 'Mt', 'autocorrthreshold', 'autocorrlag');
    n = size(data,1);
    cell_num = size(data,2);
    disp(['n = ' num2str(n) ', cell_num = ' num2str(cell_num) ', Mt = ' num2str(Mt)])
    
    Tc_real = 1:time_step:n;
    data2 = data(Tc_real, targetcells);
    
    
    %%%  clustering  %%%
    
    figure;
    Y = pdist(data2', @distfun);
    Z = linkage(Y, method);
    classnum_or_cutoff = size(data2,2);
    [H,T,outperm] = dendrogram(Z, classnum_or_cutoff); %classnumのブランチに分ける
    close(gcf)
    
    
    m = size(data2,2);
    clusters = cell(0);
    clustmat = logical(zeros(2*m-1, m));
    for i = 1:m
        clusters{i} = [i];
        clustmat(i,i) = 1;
    end
    for i = 1:m-1
        entry = [];
        clustvec = logical(zeros(1,m));
        for j = 1:2
            item = Z(i,j);
            if item<=m
                entry = [entry, item];
                clustvec(item) = true;
            else
                entry = [entry, clusters{item}];
                clustvec = clustvec | clustmat(item,:);
            end
        end
        clusters{m+i} = entry;
        if sum(clustvec) == 0
            disp(clustvec)
        end
        clustmat(m+i, :) = clustvec;
    end
    
    %%%%  find a pair of clusters, min1 and min2, that are most anticorrelated
    
    cor = corrcoef(data2);
    min1 = -1;
    min2 = -1;
    mincor = 1;
    for i = 1:2*m-2
        clustveci = clustmat(i,:);
        for j = i+1:2*m-1
            clustvecj = clustmat(j,:);
            if ~any(clustveci & clustvecj) % if no overlap in two clusters
                subcor1 = cor(clustveci, :);
                subcor2 = subcor1(:, clustvecj);
                query = sum(sum(subcor2));
                if query < mincor
                    min1 = i;
                    min2 = j;
                    mincor = query;
                end
            end
        end
    end
    %disp('end')
    %disp([min1,min2,mincor])
    %disp(clusters{min1})
    %disp(clusters{min2})
    %disp(targetcellnames(clusters{min1})')
    %disp(targetcellnames(clusters{min2})')
    
    % plot correlation matrix
    %exptitle = 'clustering_';
    %cormatrix_plot(cor(outperm,outperm), ['sample' num2str(sampleID) ' real'], 'cellname', exptitle, m, targetcellnames(outperm));%,sampleID, K, lambda);
    %cormatrix_plot(cor(outperm,outperm), ['sample' num2str(sampleID) ' real'], 'cellnumber', exptitle, m, arrayfun(@(x) num2str(x), outperm, 'UniformOutput', false));%,sampleID, K, lambda);
    
    % min1 and min2 alternate for each sample
    ABclusters{sampleID*2-1} = targetcellnames(clusters{min1});
    ABclusters{sampleID*2} = targetcellnames(clusters{min2});  
    ABclustermat(sampleID*2-1,:) = arrayfun(@(x) any(strcmp(x, targetcellnames(clusters{min1}))), conNames);
    ABclustermat(sampleID*2,:) = arrayfun(@(x) any(strcmp(x, targetcellnames(clusters{min2}))), conNames);

end

save(fullfile(common_cell_order_folder, 'ABclusters.mat'), 'ABclusters','ABclustermat');


%%% calculate overlap between selected clusters, two from each sample %%%%

disp('<< Determine A/B groups >>')
ABdist = NaN(nsamples*2);
for i = 1:2*nsamples
     ABdist(i,i) = 0;
end

for i = 1:2*nsamples-1
    for j = i+1:2*nsamples
        % count number of overlap between ith and jth member of AB
        count = sum(ABclustermat(i,:) & ABclustermat(j,:));
        dist = 1/(count+0.1);
        ABdist(i,j) = dist;
        ABdist(j,i) = dist;
    end
end


%%%%  clustering across samples 
%%%%  each member is anticorrelated clusters from each sample

figure;
Y = squareform(ABdist);
Z = linkage(Y, 'average');

[H,T,outperm] = dendrogram(Z, nsamples*2); %classnumのブランチに分ける

for nclust = 2:2*nsamples   % gradually increase cluster numbers
    T = cluster(Z,'MaxClust',nclust);
    [B I]=sort(histcounts(T,nclust), 'descend');
    disp([B(1)/B(2) B])
    if B(1)/B(2)<1.5
        break
    end
end


Aclusters = find(T==I(1));
Bclusters = find(T==I(2));
Anumber = 1;
Bnumber = 2;
if sum(ABclustermat(Aclusters, :),'all') < sum(ABclustermat(Bclusters, :),'all')
    Anumber = 2;
    Bnumber = 1;
end


%%% membershipTable: conNames x samples, value is 1 or 2.

membershipTable = zeros(length(conNames), nsamples);
for samplei = 1:nsamples
    for j=1:2
        row = samplei*2-2+j;
        if any(Aclusters == row)
            membershipTable(ABclustermat(row,:),samplei) = Anumber;
        elseif any(Bclusters == row)
            membershipTable(ABclustermat(row,:),samplei) = Bnumber;
        end
    end
end

%writematrix(membershipTable,'membershipTable.csv')

ABmembers = find(sum(membershipTable~=0, 2)>0);
numberOfOne = sum(membershipTable(ABmembers,:)==1,2);
numberOfTwo = sum(membershipTable(ABmembers,:)==2,2);
Abias = numberOfOne - numberOfTwo;
[Asort,Aorder] = sort(Abias, 'descend');


%disp(conNames(ABmembers(Aorder))')


membershipMat = [numberOfOne(Aorder) numberOfTwo(Aorder) membershipTable(ABmembers(Aorder), :)];
membershipTable2 = array2table(membershipMat);
membershipTable2.Properties.VariableNames = [{'number of ones'}, {'number of twos'}, arrayfun(@(x) ['sample' num2str(x)], 1:nsamples, 'UniformOutput', false )];
membershipTable2.Properties.RowNames = conNames(ABmembers(Aorder));

writetable(membershipTable2,[common_cell_order_folder '/membershipTable.csv'],'WriteRowNames',true)

%%% makes standard order of neurons  %%%

%ABnonmembers = 1:length(conNames);
%ABnonmembers(ABmembers) = [];
%conNamesOrder = [conNames(ABmembers(Aorder));conNames(ABnonmembers)];
conNamesOrder = conNames(ABmembers(Aorder));
save(fullfile(common_cell_order_folder, 'conNamesOrder.mat'), 'conNamesOrder')



%%% remake outperm for each sample %%%
disp('<< remake outperm >>')
samples = 1:24;
for sampleID = samples
    
    disp(['< sample ' num2str(sampleID) ' >'])
    
    load(fullfile(common_data_folder,['sample' num2str(sampleID) '_data.mat']), 'data', 'uniqNames', 'targetcellnames', 'targetcells', 'Mt', 'autocorrthreshold', 'autocorrlag');  % 
    n = size(data,1);
    cell_num = size(data2,2);
    disp(['n = ' num2str(n) ', cell_num = ' num2str(cell_num) ', Mt = ' num2str(Mt)])
    
    Tc_real = 1:time_step:n;
    data2 = data(Tc_real, targetcells);
    
    
    %%%  clustering  %%%
    
    figure;
    Y = pdist(data2', @distfun);
    Z = linkage(Y, method);
    %classnum_or_cutoff = size(data2,2);
    [H,T,outperm] = dendrogram(Z, classnum_or_cutoff); %classnumのブランチに分ける
    
    % adjust to order
    conCellOrder = [];
    for coni = 1:length(conNamesOrder)
        ui = find(strcmp(targetcellnames, conNamesOrder{coni}));
        if ~isempty(ui)
            conCellOrder = [conCellOrder, ui];
        end
    end
    targetcellrank = NaN(length(conCellOrder),1);
    targetcellrank(conCellOrder) = 1:length(conCellOrder);
    
    
    Mz = size(Z,1);
    assert(Mz == Mt-1);
    cellmembers = cell(Mt+Mz, 1);
    for i = 1:Mt
        cellmembers{i} = i;
    end
    items = NaN(1,2);cells = cell(1,2);means = NaN(1,2);
    for i = 1:Mz
        for j=1:2
            items(j) = Z(i,j);
            cells{j} = cellmembers{items(j)};
            means(j) = nanmean(-exp(-targetcellrank(cells{j})/3));
            if isnan(means(j))
                means(j) = 1;
            end
        end
        if means(2)<means(1)
            cellmembers{Mt+i} = [cells{2} cells{1}];
        else
            cellmembers{Mt+i} = [cells{1} cells{2}];
        end

    end
    outperm = cellmembers{Mt+Mz};
    
    save(fullfile(common_cell_order_folder, ['sample' num2str(sampleID) '_common_outperm.mat']), 'outperm')
    
    % plot correlation matrix
    figure;
    %exptitle = ['rearranged']; % ['sample' num2str(sampleID)]
    corrmat = corrcoef(data2);
    save_title = fullfile(common_cell_order_folder, ['sample' num2str(sampleID)]);
    cormatrix_plot(corrmat(outperm,outperm), ['sample' num2str(sampleID) ' real'], '', save_title, Mt, targetcellnames(outperm));%,sampleID, K, lambda);
    
    
end
