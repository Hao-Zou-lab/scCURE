%%%%%%%%%%% Optimize the number of clusters by using leave one out
%%%%%%%%%%% validation to maximize the prediction capability of the model

clear;
close all;

addpath('/scCURE/matlab') %set the local work path

%% load cell IDs,  patient labels and response patterns
table_pre_R = readtable('.\CD8 pre resp patientID.txt',  'Delimiter','\t', 'ReadVariableNames', 0); % from ICB database
table_pre_NR = readtable('.\CD8 pre nonresp patientID.txt',  'Delimiter','\t', 'ReadVariableNames', 0);
table_post_R = readtable('.\CD8 post resp patientID.txt',  'Delimiter','\t', 'ReadVariableNames', 0);
table_post_NR = readtable('.\CD8 post nonresp patientID.txt',  'Delimiter','\t', 'ReadVariableNames', 0);

patients_pre = {'P1','P12','P15','P2','P20','P24','P25','P26','P27','P28','P29','P3','P31','P33','P35','P4','P6','P8'};
resp_pattern = {'R', 'NR', 'NR','NR', 'NR', 'R', 'NR',  'R',  'NR',  'R', 'R',  'NR', 'NR', 'R',  'R',  'NR','NR', 'R'};

%% load scRNA-seq data
load data_CD8; %
cell_IDs_pre = strrep(cell_IDs_pre, '"', '');
cell_IDs_R_post = strrep(cell_IDs_R_post, '"', '');
cell_IDs_NR_post = strrep(cell_IDs_NR_post, '"', '');

%% leave one out validation 
AUC_all = zeros(5,5);
P_all = ones(5,5);
for m = 2:6 % the candidate cluseter numbers in group A (pre)
    for n = 2:6 % the candidate cluseter numbers in group B (post)

        itr = length(patients_pre); % LOO iterations
        nCell_R_similar = zeros(1, itr); %
        nCell_NR_similar = zeros(1, itr);
        
        for i = 1:itr % LOO for each patient
            idx = find(ismember(table2cell(table_post_R(:,2)), patients_pre{i})); % remove the test person from R_post data
            buffer_R_post = data_R_post;
            buffer_R_post(:,idx) = [];
            
            idx = find(ismember(table2cell(table_post_NR(:,2)), patients_pre{i})); % remove the test person from NR_post data
            buffer_NR_post = data_NR_post;
            buffer_NR_post(:,idx) = [];
            
            idx = find(ismember(table2cell(table_pre_R(:,2)), patients_pre{i})); % identify the pre cells from the test person
            cellIDs_pre_R = table2cell(table_pre_R(idx,1));
            idx = find(ismember(table2cell(table_pre_NR(:,2)), patients_pre{i}));
            cellIDs_pre_NR = table2cell(table_pre_NR(idx,1));
            cellIDs_pre_test = [cellIDs_pre_R, cellIDs_pre_NR];
            
            [~, ia] = intersect(cell_IDs_pre, cellIDs_pre_test); % allocate the test person cells in  pre data
            
            %%%%% R-like identification
            OGFSC_idx = OGFSC([data_pre, buffer_R_post], 'plot_option', 0, 'nBins', 20);
            data11 = data_pre(OGFSC_idx,:);
            data22 = buffer_R_post(OGFSC_idx,:);
            data = [data11, data22];
            [coeff,score] = pca(data');
            pca_score = score(:,1:5);
            data1_pca = pca_score(1:size(data11,2),:);
            data2_pca = pca_score(size(data11,2)+1:end,:);
            [Puri1, Puri2] = scCURE(data1_pca, m, data2_pca, n);
            
            nCell_R_similar(i) = length(intersect(ia, Puri1)); % calculate the number of R-like cells in the test patient
            
            %%%%% NR-like identification
            OGFSC_idx = OGFSC([data_pre, buffer_NR_post], 'plot_option', 0, 'nBins', 20);
            data11 = data_pre(OGFSC_idx,:);
            data22 = buffer_NR_post(OGFSC_idx,:);
            data = [data11, data22];
            [coeff,score] = pca(data');
            pca_score = score(:,1:5);
            data1_pca = pca_score(1:size(data11,2),:);
            data2_pca = pca_score(size(data11,2)+1:end,:);
            [Puri1, Puri2] = scCURE(data1_pca, m, data2_pca, n);
            
            nCell_NR_similar(i) = length(intersect(ia, Puri1)); % calculate the number of NR-like cells in the test patient
            
        end
        
        % save scCURE_CD8_pred nCell_NR_similar nCell_R_similar;
        
        ratio = nCell_R_similar./nCell_NR_similar; % the prediction score of each patient is the ratio of R-like/NR-like
        
        [X,Y,T,AUC] = perfcurve(resp_pattern, ratio, 'R', 'nboot', 10);% to evaluate different predictors, change line 225
        figure;
        plot(X(:,1),Y(:,1),'b','LineWidth',2);
        hold on;
        plot([0,1], [0,1], '--k');
        xlabel('1-Specificity');
        ylabel('Sensitivity');
        AUC_true = AUC(1);
        
        AUC_perm = zeros(100,1);
        for i = 1:100
            rand_idx= randperm(length(ratio));
            [X,Y,T,AUC] = perfcurve(resp_pattern(rand_idx), ratio, 'R', 'nboot', 100);
            AUC_perm(i) = AUC(1);
        end
        p = (1-cdf('norm', AUC_true, mean(AUC_perm), std(AUC_perm)));
        
        AUC_all(m-1, n-1) = AUC_true;
        P_all(m-1, n-1) = p;
    end
end

save AUC_all_CD8_LOO AUC_all P_all;

