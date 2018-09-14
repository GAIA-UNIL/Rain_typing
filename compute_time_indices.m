function[advection_Eastward, advection_Northward, corr_temporal_lag1]=compute_time_indices(matrix_data,matrix_data_tp1,step_sec,ds)
    %advection + correlation lag 1
    Im_ref=matrix_data;
    Template_to_search=matrix_data_tp1(ds:end-ds,ds:end-ds);
    C2 = normxcorr2(Template_to_search, Im_ref);
    %-----
    [sx,sy]=size(C2);
    C2=C2(floor(sx/2)-ds:floor(sx/2)+ds,floor(sy/2)-ds:floor(sy/2)+ds);
    [sx,sy]=size(C2);
    
    %-----
    [lpeak, cpeak] = find(C2==max(C2(:)));
    lpeak=mean(lpeak);
    cpeak=mean(cpeak);
    
    time_lag=1;
    dep_l = floor(lpeak-(sx/2)/time_lag);
    dep_c = floor(cpeak-(sy/2)/time_lag);
    advection_velocity=sqrt(dep_l^2+dep_c^2)*1000/(time_lag*step_sec); %!!! change the pixel resolution if needed
    
    advection_direction=atan2(dep_l,-dep_c)*180/pi; %!!!!! reference frame!!!!!!!
    sin_advection_direction=sin(advection_direction*pi/180);
    cos_advection_direction=cos(advection_direction*pi/180);
    
    advection_Eastward=advection_velocity*cos_advection_direction;
    advection_Northward=advection_velocity*sin_advection_direction;
    
    corr_temporal_lag1=max(C2(:));

end