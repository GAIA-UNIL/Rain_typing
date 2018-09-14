function[fraction_main_cluster, connectivity_index, shape_index, area_index]=compute_space_indices(grid_preprocess)
%Rq grid_preprocess is a binary matrix
    
    if (sum(sum(grid_preprocess)))<2
        fraction_main_cluster=NaN;
        connectivity_index=NaN;
        shape_index=NaN;
        area_index=NaN;
        return;
    end

    [sx,sy]=size(grid_preprocess);
    
    Mamat=(grid_preprocess>0);

    stats = regionprops(Mamat,'Area','Perimeter','ConvexArea');
    stats_field=regionprops(uint8(Mamat),'EulerNumber','ConvexArea');

    % number of rain clusters
    nb_clusters=length(cat(1,stats.Area));
    
    %fraction coverage main cluster
    fraction_main_cluster=max(cat(1,stats.Area))/(sx*sy);
    
    %connectivity index
    NP=sum(cat(1,stats.Area));
    NC=nb_clusters;
    connectivity_index=1-(NC-1)/(sqrt(NP)+NC);

    %shape index
    Perimeter_tot=sum(cat(1,stats.Perimeter));
    nb_pix_rain=sum(cat(1,stats.Area));
    Pmin=0;
    if floor(sqrt(nb_pix_rain))==sqrt(nb_pix_rain)
        Pmin=4*sqrt(nb_pix_rain);
    else
        Pmin=2*(floor(2*sqrt(nb_pix_rain))+1);
    end
    shape_index=Pmin/Perimeter_tot;

    %area index
    Area=sum(cat(1,stats.Area));
    Area_convex=stats_field.ConvexArea;
    area_index=Area/Area_convex;
    
end