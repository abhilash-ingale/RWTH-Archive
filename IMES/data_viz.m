%% Clear workspace, command window
clear;
close all;
clc;

%%
classf_feats = {'mcu_timestamp','bmi_accel_x', 'bmi_accel_y','bmi_accel_z','bmi_gyro_x','bmi_gyro_y','bmi_gyro_z'};

empty_ds = read_and_clean("empty_imu.csv",classf_feats);
wood_ds = read_and_clean("wood_imu.csv",classf_feats);
plastic_ds = read_and_clean("plastic_imu.csv",classf_feats);


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

