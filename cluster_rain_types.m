function[GMModel,Xraw_standardized,V_mean,V_std,idx]=cluster_rain_types(Xraw,nb_clusters)

[sx,sy]=size(Xraw);
Xraw_standardized=Xraw;
V_mean=[];
V_std=[];
for c=1:sy
    my_mean=mean(Xraw(:,c));
    my_std=std(Xraw(:,c));
    V_mean=[V_mean; my_mean];
    V_std=[V_std; my_std];
    if my_std>1e-5
        Xraw_standardized(:,c)=(Xraw(:,c)-my_mean)/my_std;
    else
        Xraw_standardized(:,c)=randn(sx,1)*0.001+Xraw(:,c)-my_mean;
    end
end
options = statset('MaxIter',2000,'TolFun',1e-8);
GMModel = fitgmdist(Xraw_standardized,nb_clusters,'Options',options,'CovarianceType','full','RegularizationValue',0.01);

[idx,~,~] = cluster(GMModel,Xraw_standardized);