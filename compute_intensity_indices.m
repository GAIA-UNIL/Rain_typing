function[prop_rain, mean_intensity, Q80]=compute_intensity_indices(matrix_data)
    matrix_T0=(matrix_data>0.0001);
    F=find(matrix_data>0.0001);
    V_nonzerorain=matrix_data(F);
    [sx,sy]=size(matrix_T0);
    prop_rain=sum(sum(matrix_T0))/(sx*sy);
    
    mean_intensity=mean(V_nonzerorain);
    
    Q80=quantile(V_nonzerorain,0.8);
end