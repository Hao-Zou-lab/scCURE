function [DE_up_idx, DE_down_idx] = DEG(DECells1, DECells2, nonDECells1, nonDECells2, k)


%% negative control from nonDE cells
N = size(DECells1,1);
P_nonDE = ones(1,N);
FC_nonDE = zeros(1,N);
for i = 1:N
    P_nonDE(i) = ranksum(nonDECells1(i,:), nonDECells2(i,:));
    FC_nonDE(i) = mean(nonDECells2(i,:))- mean(nonDECells1(i,:));
end
P_nonDE(find(isnan(P_nonDE))) = 1;

% P CI
P_nonDE = -log10(P_nonDE);
P_nonDE = [-P_nonDE, P_nonDE]; % Mean Of All Experiments At Each Value Of ‘x’
P_nonDECI = k*std(P_nonDE);


% positive FC CI
FCp = FC_nonDE(FC_nonDE>0);
FCp = [-FCp, FCp];
FCp_nonDECI = k*std(FCp);

% negative FC CI
FCn = FC_nonDE(FC_nonDE<0);
FCn = [-FCn, FCn];
FCn_nonDECI = k*std(FCn);

%% DE gene identification from DE cells

P = ones(1,N);
FC = zeros(1,N);
for i = 1:N
    P(i) = ranksum(DECells1(i,:), DECells2(i,:));
    FC(i) = mean(DECells2(i,:))- mean(DECells1(i,:));
end
P(find(isnan(P))) = 1;
P = -log10(P);
idx_p = find(P>P_nonDECI);

DE_up_idx = intersect(idx_p, find(FC>FCp_nonDECI));
DE_down_idx = intersect(idx_p, find(FC<-FCn_nonDECI));
