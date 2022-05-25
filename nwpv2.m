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
d=1-normcdf(auto(-norminv(d)));
%d=1-normcdf(-norminv(d));
%combining the individual p values.
[y,z]=realcomp(d,opt);
if isempty(alpha)
    alpha=0.05;
end
ind = find(y<alpha);
if ls==1
    vidz(d,opt,ind);
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
