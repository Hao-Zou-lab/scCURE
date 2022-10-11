function [Unchanged1, Unchanged2] = scCURE(data1, K1, data2, K2)

%%%%  Input data1, data2 - dimentional reduced data matrix with PCs vs cells from two groups.
%%%%        K1, K2 - the numbers of clusters 
%%%%  Output Unchanged1, Unchanged2 - the steady cells from two groups

idx1 = [];
idx2 = [];
II = 20; % repeat for 20 iterations
for I = 1:II
    %% build GMM
    GMModel_1 = fitgmdist(data1, K1, 'CovarianceType','full', 'Replicates',5, 'RegularizationValue',0.00001);
    GMModel_2 = fitgmdist(data2, K2, 'CovarianceType','full', 'Replicates',5, 'RegularizationValue',0.00001);
    Sigma1 = squeeze(GMModel_1.Sigma);
    Sigma2 = squeeze(GMModel_2.Sigma);
    
    %% calculate mutual Kullback-Leibler divergence
    D12 = zeros(K1,K2);
    for i = 1:K1
        mu1 = GMModel_1.mu(i,:);
        sig1 = Sigma1(:,:,i);
        for j = 1:K2
            mu2 = GMModel_2.mu(j,:);
            sig2 = Sigma2(:,:,j);
            D12(i,j) = mvgkl(mu1', sig1, mu2', sig2);
        end
    end
    
    %% find the mutually cloest Gaussian models
    [~, IX] = min(D12,[],2);
    [~, IX2] = min(D12,[],1);
    Pair12 = [];
    for i = 1:K1
        C = IX(i);
        buffer = IX2(C);
        if i==buffer
            Pair12 = [Pair12; [i,C]];
        end
    end
    
    %% assign cells to specific Gaussian models
    cluster1 = cluster(GMModel_1,data1); % Cluster index
    cluster2 = cluster(GMModel_2,data2); % Cluster index
    
    %% identify the unchanged cells based on the mutually cloest clustres
    idx1_temp = [];
    idx2_temp = [];
    for i = 1:size(Pair12,1)
        idx1_temp = [idx1_temp; find(cluster1==Pair12(i,1))];
        idx2_temp = [idx2_temp; find(cluster2==Pair12(i,2))];
    end
    
    idx1 = [idx1; idx1_temp];
    idx2 = [idx2; idx2_temp];
end

%% the cells identified as unchanged for more than 10 out of 20 times are finally identified as unchanged
T1 = tabulate(idx1);
idx = find(T1(:,2)>10);
Unchanged1 = T1(idx,1);

T2 = tabulate(idx2);
idx = find(T2(:,2)>10);
Unchanged2 = T2(idx,1);
