function Result=DEGstat_20170415(rT,rF,cont,treat,n_treat,t_opt)
rT(find(rT==inf))=max(max(rT(rT~=inf),[],1),2);  
rT(rT==-inf)=min(min(rT(rT~=-inf),[],1),2);
rT(isnan(rT)==1)=0;

% t_opt=1 : one sample and paired t-test
% t_opt=2 : two sample t-test, unequal variance

%% Calculation of p-values for t-test and fold change test
if t_opt==1    
    for j=1:(size(treat,2)/n_treat)
        [~, ~, ~, t1]=ttest(treat(:,(n_treat*(j-1)+1):(n_treat*j)),cont,[],'both',2);
        Result.T(:,j)=t1.tstat;
        Result.T(isnan(Result.T(:,j))==1,j)=0;      % Expression이 전부 0인경우 T-test결과 NaN이 나옴
        Result.T(find(Result.T(:,j)==inf),j)=max(Result.T(Result.T(:,j)~=inf,j),[],1);       % T=inf의 경우 max(T)로
        Result.T(find(Result.T(:,j)==-inf),j)=min(Result.T(Result.T(:,j)~=-inf,j),[],1);      % T=-inf의 경우 min(T)로
        Result.Pt(:,j)=pval2tail(rT(:),Result.T(:,j));
        
        Result.FC(:,j)=nanmedian(treat(:,(n_treat*(j-1)+1):(n_treat*j))-cont,2);
        Result.Pf(:,j)=pval2tail(rF(:),Result.FC(:,j));
        
        clear t1;
    end

elseif t_opt==2    
    for j=1:(size(treat,2)/n_treat)
        [~, ~, ~, t2]=ttest2(treat(:,(n_treat*(j-1)+1):(n_treat*j)),cont,[],'both','unequal',2);
        Result.T(:,j)=t2.tstat;
        Result.T(isnan(Result.T(:,j))==1,j)=0;      % Expression이 전부 0인경우 T-test결과 NaN이 나옴
        Result.T(find(Result.T(:,j)==inf),j)=max(Result.T(Result.T(:,j)~=inf,j),[],1);       % T=inf의 경우 max(T)로
        Result.T(find(Result.T(:,j)==-inf),j)=min(Result.T(Result.T(:,j)~=-inf,j),[],1);      % T=-inf의 경우 min(T)로
        Result.Pt(:,j)=pval2tail(rT(:),Result.T(:,j));
        
        Result.FC(:,j)=nanmedian(treat(:,(n_treat*(j-1)+1):(n_treat*j)),2)-nanmedian(cont,2);
        Result.Pf(:,j)=pval2tail(rF(:),Result.FC(:,j));
        clear t2;
    end
end

%% Integration of P-values
for j=1:(size(treat,2)/n_treat)
    d(:,1) = 1-normcdf(auto(-norminv(Result.Pt(:,j))));
    d(:,2) = 1-normcdf(auto(-norminv(Result.Pf(:,j))));
    [pcom,~,~]=nwpv2(d,3,0.05,0);
    Result.Pcom(:,j)=pcom;
    clear d;
end


Result.cont=cont;
Result.treat=treat;
Result.t_opt=t_opt;
end


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

end



%-------------------------------------------------------------------------------
% Copyright (C) 2005 by Institute for Systems Biology,
% Seattle, Washington, USA.  All rights reserved.
% 
% This source code is distributed under the GNU Lesser 
% General Public License, the text of which is available at:
%   http://www.gnu.org/copyleft/lesser.html
%
%-------------------------------------------------------------------------------
% Module:    nwpv2.m
%
% Author:    Daehee Hwang 
%            Institute for Systems Biology
%
% Function:  non-weighted integration methods
%-------------------------------------------------------------------------------

function [y,ind,z]=nwpv2(d,opt,alpha,ls);

%correction to ensure numerical instability
[n,k]=size(d);
for i=1:k
    tn=find(d(:,i)==0);
    if ~isempty(tn)
        d(find(d(:,i)==0),i)=linspace(1.e-10,1.e-323,length(tn))';
    end
    if ~isempty(find(d(:,i)<1.e-323))
        d(find(d(:,i)<1.e-323),i)=min(d(find(d(:,i)>1.e-323),i))
    end
    if ~isempty(find(d(:,i)>0.99999999999999994))
        d(find(d(:,i)>0.99999999999999994),i)=max(d(find(d(:,i)<0.99999999999999994),i));
    end
end
%normalization to treat the individual vaiables equally.
%d=1-normcdf(auto(-norminv(d)));
d=1-normcdf(-norminv(d));
%combining the individual p values.
[y,z]=realcomp(d,opt);
if isempty(alpha)
    alpha=0.05;
end
ind = find(y<alpha);
if ls==1
    vidz(d,opt,ind);
end
end

function [y,z]=realcomp(x,opt);
[n,k]=size(x);
%multiplication const for Fisher's, MG, and Stouffer's methods.
if opt==1
    z=-2*log(x);
    k=sum(~isnan(z),2);
    z=nansum(z,2);
    y=1-chi2cdf(z,2*k);
elseif opt==2
    z1=log(x./(1-x));
    k=sum(~isnan(z1),2);
    ind=find(k>0);
    z=NaN*ones(size(z1,1),1);
    z(ind)=nansum(z1(ind,:),2).*(-sqrt((15*k(ind)+12)./((5*k(ind)+2).*k(ind)*pi^2)));
    y=1-tcdf(z,5*k+4);
elseif opt==3
    z1=-norminv(x);z=NaN*ones(size(z1,1),1);
    k=sum(~isnan(z1),2);
    ind=find(k>0);
    z=NaN*ones(size(z1,1),1);
    z(ind)=nansum(z1(ind,:),2)./sqrt(k(ind));
    y=1-normcdf(z);    
end
end

function vidz(d,opt,ind);
%pseudo P variable generation
if size(d,2)>2
    px=realcomp(d(:,1:end-1),opt);
else
    px=d(:,1);    
end
plot(px,d(:,end),'g.'); hold on;
plot(px(ind),d(ind,end),'mo');
axis([0 1 0 1]);
end