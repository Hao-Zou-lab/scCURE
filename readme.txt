clear;

close all;





addpath('/scCURE/matlab') %set the local work path

%% 



% fileDir = '.\CD8 pre resp.xls'; %read entire file into string

"% [data_R_pre, geneList, cell_types1, cell_IDs_unique1, cell_IDs_R_pre] = readData(fileDir);"

% 

% fileDir = '.\CD8 pre nonresp.xls'; %read entire file into string

"% [data_NR_pre, geneList, cell_types1, cell_IDs_unique1, cell_IDs_NR_pre] = readData(fileDir);"

% 

"% data_pre = [data_R_pre, data_NR_pre];"

"% cell_IDs_pre = [cell_IDs_R_pre, cell_IDs_NR_pre];"

% 

% fileDir = '.\CD8 post resp.xls'; %read entire file into string

"% [data_R_post, geneList, cell_types1, cell_IDs_unique1, cell_IDs_R_post] = readData(fileDir);"

% 

% fileDir = '.\CD8 post nonresp.xls'; %read entire file into string

"% [data_NR_post, geneList, cell_types1, cell_IDs_unique1, cell_IDs_NR_post] = readData(fileDir);"

% % 

% save data_CD8 data_pre cell_IDs_pre data_R_post cell_IDs_R_post data_NR_post cell_IDs_NR_post geneList;

%% Responders



%%%%%%%%% load scRNA-seq data

load data_CD8;



%%%%%%%%% gene filtering for random noise reduction using OGFSC

"OGFSC_idx = OGFSC([data_pre, data_R_post], 'plot_option', 1, 'nBins', 20);"

"data11 = data_pre(OGFSC_idx,:); % pre-treatment sample"

"data22 = data_R_post(OGFSC_idx,:); % post-treatment"



%%%%%%%%% PCA for dimentional reduction

"data = [data11, data22];"

"[coeff,score] = pca(data);"

"pca_score = score(:,1:5);"

"data1_pca = pca_score(1:size(data11,2),:);"

"data2_pca = pca_score(size(data11,2)+1:end,:);"



%%%%%%%%% (optional) estimate the number of Gaussian models by AIC

estimateK(data1_pca);

estimateK(data2_pca);



%%%%%%%%% scCURE 

"[Puri1, Puri2] = scCURE(data1_pca, 5, data2_pca, 3);"



%%%%%%%%% put 'changed' or 'steady' label for each cell

"opt1 = repmat({'changed'}, 1,size(data1_pca,1));"

opt1(Puri1) = {'steady'};



"opt2 = repmat({'changed'}, 1,size(data2_pca,1));"

opt2(Puri2) = {'steady'};



%%%%%%%%%% save outputs

"cellSamples = [repmat({'Pre'},1,size(data1_pca,1)), repmat({'Post'}, 1,size(data2_pca,1))];"

geneList_filtered = geneList(OGFSC_idx);

"opt = [opt1, opt2];"

save data_processed_R data11 data22 geneList_filtered cellSamples opt data_pre data_R_post geneList cell_IDs_pre cell_IDs_R_post;









%% non-Responders

%%%%%%%%% load scRNA-seq data

load data_CD8;



%%%%%%%%% gene filtering for random noise reduction using OGFSC

"OGFSC_idx = OGFSC([data_pre, data_NR_post], 'plot_option', 1, 'nBins', 20);"

"data11 = data_pre(OGFSC_idx,:);"

"data22 = data_NR_post(OGFSC_idx,:);"



%%%%%%%%% PCA for dimentional reduction

"data = [data11, data22];"

"[coeff,score] = pca(data');"

"pca_score = score(:,1:5);"

"data1_pca = pca_score(1:size(data11,2),:);"

"data2_pca = pca_score(size(data11,2)+1:end,:);"



%%%%%%%%% (optional) estimate the number of Gaussian models by AIC

estimateK(data1_pca);

estimateK(data2_pca);



%%%%%%%%% scCURE 

"[Puri1, Puri2] = scCURE(data1_pca, 5, data2_pca, 3);"





%%%%%%%%% put 'changed' or 'steady' label for each cell

"opt1 = repmat({'changed'}, 1,size(data1_pca,1));"

opt1(Puri1) = {'steady'};



"opt2 = repmat({'changed'}, 1,size(data2_pca,1));"

opt2(Puri2) = {'steady'};





%%%%%%%%%% save outputs

"cellSamples = [repmat({'Pre'},1,size(data1_pca,1)), repmat({'Post'}, 1,size(data2_pca,1))];"

geneList_filtered = geneList(OGFSC_idx);

"opt = [opt1, opt2];"

save data_processed_NR data11 data22 geneList_filtered cellSamples opt data_pre data_NR_post geneList cell_IDs_pre cell_IDs_NR_post;