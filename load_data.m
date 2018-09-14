function[rain_data_matrix]=load_data(sec)
    %The radar images must be croped beforehand to cover only the area of interest
    %The radar images must be in mm/h
    cd('Data\')    
    my_filename=strcat('Data_synthetic_radarserie_',num2str(floor(sec/(10*60))),'.mat'); %Here we consider one radar image every 10min
    rain_data_matrix=load(my_filename);
    rain_data_matrix=rain_data_matrix.Mamat;
    cd('..')
end