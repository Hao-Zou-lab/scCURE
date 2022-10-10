function estimateK(data)


K_potential = 2:16;
GMModel = cell(1,length(K_potential));
aic = zeros(1,length(K_potential));
bic = zeros(1,length(K_potential));

for k = 1:length(K_potential)
    GMModel{k} = fitgmdist(data, K_potential(k), 'CovarianceType','full', 'RegularizationValue',0.001);
    aic(k) = GMModel{k}.AIC;
%     bic(k) = GMModel{k}.BIC;
end

figure;
bar((aic));
ylabel('AIC');
xlabel('numComponents');
xticks(1:length(K_potential));
xticklabels(K_potential);
% figure;
% bar((bic));
% ylabel('BIC');
