function[V_raintypes_segment]=segment_rain_events(idx,V_corresp_idx,V_time,V_time_all,length_min_segment,dt_inter_events)

%split into rain events
ind_struct=1;
ind_start=1;

for i=2:length(V_time)
    t1=V_time(i-1);
    t2=V_time(i);
    dt = t2-t1;
    if dt>dt_inter_events
        Strct_rain_event(ind_struct).time=V_time(ind_start:i-1);
        ind_start=i;
        ind_struct=ind_struct+1;
    end
end
Strct_rain_event(ind_struct).time=V_time(ind_start:i);

%apply GMM classification of relevant (enough rain in corresponding radar image)
for i=1:length(Strct_rain_event)
    for j=1:length(Strct_rain_event(i).time)
        Strct_rain_event(i).raw_rain_type(j)=0;
        my_ind=find(V_time==Strct_rain_event(i).time(j));
        if isempty(my_ind)==0
            if V_corresp_idx(my_ind)>0
                Strct_rain_event(i).raw_rain_type(j)=idx(my_ind);
            end
        end
    end
end


%filter noise
se1 = strel('line',length_min_segment,0);
my_nb_clusters=max(idx);
for i=1:length(Strct_rain_event)
    Strct_rain_event(i).rain_type_filtred=Strct_rain_event(i).raw_rain_type*0.0;
    for j=1:my_nb_clusters
        MMM=(Strct_rain_event(i).raw_rain_type==j);
        MMM=imopen(MMM,se1);
        for k=1:length(MMM)
            if MMM(k)>0
                Strct_rain_event(i).rain_type_filtred(k)=MMM(k)*j;
            end
        end
    end
end

%fill zeros with nearest neighboor
for i=1:length(Strct_rain_event)
    if sum(Strct_rain_event(i).rain_type_filtred)>0
        %add first / last value to avoid NaN in nearest neighbor interpolation
        V_ind_nonzero=find(Strct_rain_event(i).rain_type_filtred>0);
        Strct_rain_event(i).rain_type_filtred(1)=Strct_rain_event(i).rain_type_filtred(V_ind_nonzero(1));
        Strct_rain_event(i).rain_type_filtred(end)=Strct_rain_event(i).rain_type_filtred(V_ind_nonzero(end));
        
        %nearest neighbor interpolation
        Strct_rain_event(i).rain_type_final=Strct_rain_event(i).rain_type_filtred;
        xx=1:1:length(Strct_rain_event(i).time);
        vv=Strct_rain_event(i).rain_type_filtred;
        V_ind=find(vv>0);
        xd=xx(V_ind);
        vd=vv(V_ind);
        V_ind_xq=find(vv==0);
        xq=xx(V_ind_xq);
        vq=interp1(xd,vd,xq,'nearest');
        Strct_rain_event(i).rain_type_final(V_ind_xq)=vq;
    else
        ind_nonzero=find(Strct_rain_event(i).raw_rain_type>0);
        Strct_rain_event(i).rain_type_final=ones(1,length(Strct_rain_event(i).raw_rain_type))*mode(Strct_rain_event(i).raw_rain_type(ind_nonzero));
    end
end

%final structure
V_raintypes_segment=zeros(length(V_time_all),1);
for i=1:length(Strct_rain_event)
    for j=1:length(Strct_rain_event(i).time)
        my_ind=find(V_time_all==Strct_rain_event(i).time(j));
        if ~isempty(my_ind) && ~isnan(my_ind) && ~isnan(Strct_rain_event(i).rain_type_final(j))
            V_raintypes_segment(my_ind)=Strct_rain_event(i).rain_type_final(j);
        end
    end
end

end