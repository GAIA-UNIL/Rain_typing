function[M_BIC]=cluster_BIC(X_raw,max_clusters)


[sx,sy]=size(X_raw);
for c=1:sy
    my_mean=mean(X_raw(:,c));
    my_std=std(X_raw(:,c));
    if my_std>1e-5
        X_raw(:,c)=(X_raw(:,c)-my_mean)/my_std;
    else
        X_raw(:,c)=randn(sx,1)*0.001+X_raw(:,c)-my_mean;
    end
end

M_BIC=[];
options = statset('MaxIter',2000,'TolFun',1e-8);
for j=1:max_clusters
    GMModel = fitgmdist(X_raw,j,'Options',options,'CovarianceType','full','RegularizationValue',0.1);
    M_BIC=[M_BIC;[j GMModel.BIC GMModel.AIC GMModel.NegativeLogLikelihood]];
end

end