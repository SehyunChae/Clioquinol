function [rT rF]=emp_dist_v3(cont,treat,n_cont,n_treat,perm,opt,t_opt)
% opt=1 : random permutation
% opt=2 : random sampling

% t_opt=1 : one sample and paired t-test
% t_opt=2 : two sample t-test, unequal variance

%% Empirical distribution of log2 fold changes
log2_normdata=[cont treat];
n_cont_type=size(cont,2)/n_cont;

if t_opt==1

    if opt==1

        for i=1:perm
            idx=randperm(size(log2_normdata,2));
      
            pseudo_cont=log2_normdata(:,idx(1:n_cont));
            pseudo_treat=log2_normdata(:,idx(n_cont*n_cont_type+1:end));

            pseudo_dif=pseudo_treat-repmat(pseudo_cont,[1 size(treat,2)/n_treat]);
            medif=zeros(size(treat,1),size(treat,2)/n_treat);
            
            pseudo_t=zeros(size(treat,1),size(treat,2)/n_treat);
            
            for j=1:size(treat,2)/n_treat
                medif(:,j)=nanmedian(pseudo_dif(:,(n_treat*(j-1)+1):(n_treat*j)),2);
                [~, ~, ~, t1]=ttest(pseudo_treat(:,(n_treat*(j-1)+1):(n_treat*j)),pseudo_cont,[],'both',2);
                pseudo_t(:,j)=t1.tstat;
                clear t1
            end
                
            rT(:,i)=pseudo_t(:);
            rF(:,i)=medif(:);
            clear pseudo_cont pseudo_treat pseudo_dif medif pseudo_t;
        end

    elseif opt==2
        for i=1:perm
            idx=randperm(size(log2_normdata,2),(n_cont+n_treat));
                
            pseudo_cont=log2_normdata(:,idx(1:n_cont));
            pseudo_treat=log2_normdata(:,idx(n_cont+1:end));

            pseudo_dif=pseudo_treat-pseudo_cont;
            medif=nanmedian(pseudo_dif,2);
            [~, ~, ~, t1]=ttest(pseudo_treat,pseudo_cont,[],'both',2);
            pseudo_t=t1.tstat;
            rT(:,i)=pseudo_t;
            rF(:,i)=medif;
            clear t1 pseudo_cont pseudo_treat pseudo_dif medif pseudo_t;
        end
    
    end


elseif t_opt==2

    if opt==1

        for i=1:perm
            idx=randperm(size(log2_normdata,2));
            pseudo_med_cont(:,1)=nanmedian(log2_normdata(:,idx(1:n_cont)),2);
            pseudo_med_treat=zeros(size(treat,1),size(treat,2)/n_treat);
                
            pseudo_cont=log2_normdata(:,idx(1:n_cont));
            pseudo_treat=log2_normdata(:,idx(n_cont*n_cont_type+1:end));
            
            pseudo_t=zeros(size(treat,1),size(treat,2)/n_treat);
    
            for j=1:size(treat,2)/n_treat
                pseudo_med_treat(:,j)=nanmedian(pseudo_treat(:,(n_treat*(j-1)+1):(n_treat*j)),2);
                [~, ~, ~, t2]=ttest2(pseudo_treat(:,(n_treat*(j-1)+1):(n_treat*j)),pseudo_cont,[],'both','unequal',2);
                pseudo_t(:,j)=t2.tstat;
                clear t2
            end
            medif=pseudo_med_treat-repmat(pseudo_med_cont,[1 size(treat,2)/n_treat]);
                
            rT(:,i)=pseudo_t(:);
            rF(:,i)=medif(:);
            clear pseudo_med_cont pseudo_cont pseudo_treat pseudo_med_treat medif pseudo_t;
        end

    elseif opt==2
        for i=1:perm
            idx=randperm(size(log2_normdata,2),(n_cont+n_treat));
            pseudo_med_cont(:,1)=nanmedian(log2_normdata(:,idx(1:n_cont)),2);
            pseudo_med_treat(:,1)=nanmedian(log2_normdata(:,idx(n_cont+1:end)),2);
                
            pseudo_cont=log2_normdata(:,idx(1:n_cont));
            pseudo_treat=log2_normdata(:,idx(n_cont+1:end));

            medif=pseudo_med_treat-pseudo_med_cont;
            [~, ~, ~, t2]=ttest2(pseudo_treat,pseudo_cont,[],'both','unequal',2);
            pseudo_t=t2.tstat;
            rT(:,i)=pseudo_t;
            rF(:,i)=medif;
            clear t2 pseudo_med_cont pseudo_cont pseudo_treat pseudo_med_treat medif pseudo_t;
        end
    
    end
end
end