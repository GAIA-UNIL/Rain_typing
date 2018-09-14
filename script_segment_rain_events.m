close all
clear all
display('start processing')
%*************************Parameters of the script - to change according to the application*************************
ep_ini=1; %# of the first image to process
ep_end=856; %# of the last image to process
step_sec=10*60; %time step (in seconds) between two subsequent radar images
threshold_size=5; %minimum number of connex pixel to form a rainy patch
threshold_to_be_rainy_epoch=5; %minimum number of rainy pixels in an image to consider that the image is 'rainy'
ds=2; %search window for the correlation between two subsequent images
max_clusters=6; %maximum number of clusters tested to find the number of gaussian mixtures used in the GMM model
%Data files must be placed in the \Data folder
%The function load_data must be modified in regard to the name and structure of the data files
%*******************************************************************************************************************

%---------------------Step 1: compute indices--------------
sec_ini=ep_ini*10*60;
sec_end=ep_end*10*60;
Matrix_indices=zeros(10,ep_end-ep_ini)*NaN;
ind_column=0;

[matrix_data_tp1]=load_data(sec_ini); %The function load_data must be modified in regard to the name and structure of the data files
for sec=sec_ini:step_sec:sec_end-step_sec
    
    %load data
    ind_column=ind_column+1;
    matrix_data=matrix_data_tp1;
    [rain_data_matrix_tp1]=load_data(sec+step_sec);
    
    %pre-processing: remove connex components <threshold_size (i.e. possible groung clutters)
    M_indic_rain_tp1=(rain_data_matrix_tp1>=0.0001);
    CC = bwconncomp(M_indic_rain_tp1);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx=find(numPixels<=threshold_size);
    for k=1:length(idx)
        M_indic_rain_tp1(CC.PixelIdxList{idx(k)}) = 0;
    end
    matrix_data_tp1=rain_data_matrix_tp1.*M_indic_rain_tp1;
    if sum(sum(matrix_data>0.0001))<threshold_to_be_rainy_epoch || sum(sum(matrix_data_tp1>0.0001))<threshold_to_be_rainy_epoch
        matrix_data_tm1=matrix_data;
        matrix_data=matrix_data_tp1;
        continue;
    end
    
    %compute and store indices
    [Matrix_indices(1,ind_column),Matrix_indices(2,ind_column),Matrix_indices(3,ind_column)]=compute_intensity_indices(matrix_data);
    matrix_T0=(matrix_data>0.1);%generate binary image for Spatial Indices (SI) computation (truncation threshold = 0.1)
    [Matrix_indices(4,ind_column),Matrix_indices(5,ind_column),Matrix_indices(6,ind_column),Matrix_indices(7,ind_column)]=compute_space_indices(matrix_T0);
    [Matrix_indices(8,ind_column),Matrix_indices(9,ind_column),Matrix_indices(10,ind_column)]=compute_time_indices(matrix_data,matrix_data_tp1,step_sec,ds);
    
end
display('Step 1 [compute indices] completed')
%%
%-------------------------Step 2: Extract indices for time steps with more than 10% rainy pixels-----------------------------
prop_rain=0.1; %minimum rain fraction per time step to compute the indices
V_time=[];
V_time_all=[];
Xraw=[];
M_Results_marginal=[];
M_Results_spatial_structure=[];
M_Results_temporal_structure=[];
V_corresp_idx=ones(length(Matrix_indices(1,:)),1)*-1;
ind_corresp=1;
for i=1:length(Matrix_indices(1,:))
    V_time_all=[V_time_all;sec_ini+i*60*10];
    if sum(isinf(Matrix_indices(:,i)))==0 && ~isnan(Matrix_indices(1,i))
        if Matrix_indices(1,i)>prop_rain
            V_time=[V_time;sec_ini+i*60*10];
            Xraw=[Xraw;Matrix_indices(:,i)'];
            V_corresp_idx(i)=ind_corresp;
            ind_corresp=ind_corresp+1;
        end
        if Matrix_indices(1,i)>0.05 && Matrix_indices(1,i)<=prop_rain
            V_corresp_idx(i)=0;
        end
    end
end
display('Step 2 [extract indices] completed')

%%
%-----------------Step 3 [optional]: find the 'best' number of components using BIC criterion-------------------------
[M_BIC]=cluster_BIC(Xraw,max_clusters);
figure(1)
plot(M_BIC(:,1),M_BIC(:,2))
xlabel('number of clusters')
ylabel('BIC')
display('Step 3 [select number of gaussian mixtures] completed')

%%
%-----------------Step 4: classification using GMM-------------------------------------- 
nb_clusters=4;%should be changed according to the result of the previous step
[GMModel,Xraw_standardized,V_mean,V_std,idx]=cluster_rain_types(Xraw,nb_clusters);
%Figure 2: histograms (diagonal) and bivariate distributions of the indices. Red lines correspond to the model that has been fit
plot_index_data_standardized(Xraw_standardized,2,GMModel)
display('Step 4 [classification] completed')


%-----------------Step 5: clean the classification result and derive final rain types-----------------------
%change segment_rain_events!!!
[V_raintypes]=segment_rain_events(idx,V_corresp_idx,V_time,V_time_all,6,30*60);
%Figure 3: Rain types
figure(3)
plot(V_time_all,V_raintypes,'.')
xlabel('time (sec)')
ylabel('rain type')
display('Step 5 [Rain typing] completed')
display('End processing')
