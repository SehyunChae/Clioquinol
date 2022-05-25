% P-value computation routine
%
% p-value = pval2tail(so,s)
%
% s0 : null statistic, s : observed statistic
% p-value : computed by 2-tail test
%
% Last modified on 17th, July, 2008 by jhcho

function p = pval2tail(s0,s)

[f0,s0] = ecdf(s0); % f0 : cdf. value of null statistics, s0 : null statistics
p = interp1(s0(2:end),f0(2:end),s);
% If a certain value of s is larger/smaller than s0, interp1 returns 'NaN'
% In this case, p-value reaches 0
p(isnan(p)==1) = 0; 
p = 2*min([p 1-p],[],2);

% Numerical stabilization for p-value integration
% Removing 0 and 1 for stable norminv transformation
p(p==0) = min(p(p~=0))/2;
p(p==1) = (1-max(p(p~=1)))/2+max(p(p~=1));