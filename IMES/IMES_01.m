%% Clear variables and windows
clear;
close all;
clc;

%% READ AND PLOT DATA

classf_feats = {'mcu_timestamp','bmi_accel_x', 'bmi_accel_y','bmi_accel_z','bmi_gyro_x','bmi_gyro_y','bmi_gyro_z'};
%ds_names = ["empty_ds", "wood_ds", "plastic_ds"];

empty_ds = read_and_clean("empty_imu.csv",classf_feats);
wood_ds = read_and_clean("wood_imu.csv",classf_feats);
plastic_ds = read_and_clean("plastic_imu.csv",classf_feats);

plot_ds(empty_ds);
plot_ds(wood_ds);
plot_ds(plastic_ds);

%% APPENDIX: User-defined functions

% 1. read_and_clean : 

function temp_ds_array = read_and_clean(file_name,classf_feats)
    
    % Import .csv files into datastore
    temp_ds = datastore(file_name,'TreatAsMissing','NA', 'MissingValue',0);
    
    % Convert selected features from datastore to array 
    temp_ds.SelectedVariableNames = classf_feats;
    temp_ds_array = table2array(readall(temp_ds));
    
    % Normalize data to start from 0
    temp_ds_array(:,1) = temp_ds_array(:,1)- temp_ds_array(1,1);
    
    % Ignore rows with timestamp = 0 
    posn = find(temp_ds_array(:,1)==0);
    temp_ds_array(posn,:)= [];
    
end

% 2. plot_ds
function plot_ds(temp_ds)
    % Define labels for each sub-plot
    labels = ["bmi_accel_x", "bmi_accel_y","bmi_accel_z","bmi_gyro_x","bmi_gyro_y","bmi_gyro_z"];
    figure;
    
    % Plot accelerometer data
    for i1 = 2:4    
        j1 = i1-1 ;
        subplot(2,3,j1);
        plot(temp_ds(:,i1))
        title(labels(j1));
    end
    
    % Plot gyro data
    for i2 = 5:7    
        j2 = i2-1 ;
        subplot(2,3,j2);
        plot(temp_ds(:,i2))
        title(labels(j2));
    end
end