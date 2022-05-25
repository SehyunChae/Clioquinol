%% Data import (probe intensities from microarray)
[data, result]= readtext('Input.txt', '\t', '', '', 'numeric');

%% Quantile-normalization
data = log2(data);
norm = quantilenorm(data);

%% Null distribution generation
%[rT rF]=emp_dist_v3(Control,Treat,n_cont,n_treat,perm,opt,t_opt)
% n_cont : number of "control" samples
% n_treat : number of "treat" samples
% perm: number of permutations (e.g. 1000)
% opt=1 : random permutation
% opt=2 : random sampling
% t_opt=1 : one sample and paired t-test
% t_opt=2 : two sample t-test, unequal variance

[rT rF]=emp_dist_v3(norm(:,1:2),norm(:,3:6),2,2,1000,2,2)

%% Statistical test
% CQ vs. normoxia
Result=DEGstat_20170415(rT,rF,norm(:,1:2),norm(:,3:4),2,2)

% hypoxia vs. normoxia
Result=DEGstat_20170415(rT,rF,norm(:,1:2),norm(:,5:6),2,2)
