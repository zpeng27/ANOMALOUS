load Disney.mat;

[n,~] = size(X);
X = normalizeFea(X, 0);
Xt = X';

niters = 100;
alpha = 0.015;
beta = 0.01;
gamma = 0.009;
phi = 0.6;

At = A';
Anew = max(A,At);
L = computelaplacian(Anew, 'undirected');
R = anomalous(Xt, A, L, alpha, beta, gamma, phi,niters);
Rt = R';
score = sum(Rt.*Rt,2);
[~,idx] = sort(score, 'descend');

gnd_data = zeros(n,2);
gnd_data(:,1) = gnd;
gnd_data(:,2) = score;
[tp,fp] = roc([gnd_data(:,1),gnd_data(:,2)]);
% plot(fp,tp);
auc_value = auc(gnd_data);
fprintf('AUC_value: %f',auc_value);
% disp(score);
% disp(idx);