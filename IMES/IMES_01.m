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

% Visualize data 
% Note - X-axis shows no. of points instead of time since all points are sampled at equal frequency i.e 1600 Hz
%plot_ds(empty_ds);
%plot_ds(wood_ds);
%plot_ds(plastic_ds);

%% PREPARE DATA FOR ML

% Convert dataset array to the corresponding feature matrix
empty_ds_ml = ds_with_features(empty_ds,400,3);
wood_ds_ml = ds_with_features(empty_ds,400,3);
plastic_ds_ml = ds_with_features(empty_ds,400,3);

%% APPENDIX: User-defined functions

% 1. read_and_clean : 
% WHAT IT DOES - returns the dataset as an array reading required features from csv file
% file_name : string, name of the csv file i.e. example.csv
% classf_feats : a vector containing required feature names as strings

function temp_ds_array = read_and_clean(file_name,classf_feats)
    
    % Import .csv files into datastore
    temp_ds = datastore(file_name,'TreatAsMissing','NA', 'MissingValue',0);
    
    % Convert selected features from datastore to array 
    temp_ds.SelectedVariableNames = classf_feats;
    temp_ds_array = table2array(readall(temp_ds));
    
    % Start timestamp counter from 0
    temp_ds_array(:,1) = temp_ds_array(:,1)- temp_ds_array(1,1);
    
    % Ignore rows with timestamp = 0
    % Why? - They do not have all the 1600 points for the 1st second (observed in all 3 datasets)
    posn = temp_ds_array(:,1)==0;
    temp_ds_array(posn,:)= [];
 
end

% 2. plot_ds :
% WHAT IT DOES - returns a 6 i.e. ax_ay,az,gx,gy,gz quantity sub-plot figure

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

% 3. ds_with_features :

% WHAT IT DOES - returns a feature matrix given a dataset
% temp_ds : First argument is the dataset itself,
% sampling_factor: second argument is the no. of samples that we take for calculating features 
% model_feats: third argument is no.of features (for example, set model_feats=3 for mean,rms,amp feats)

function sampled_ds = ds_with_features(temp_ds, sampling_factor, model_feats)
    
    % Remove the timestamp column
    temp_ds = temp_ds(:,2:7);

    % Get the length of the dataset
    rows = size(temp_ds,1);
    
    % Calculate total no. of rows to be used
    sample_count = floor(rows/sampling_factor);
    
    % Remove extra rows to make no. of rows perfectly divisible by sampling factor
    temp_ds = temp_ds(1:sample_count*sampling_factor,:);
    
    % Initialize
    sampled_ds = zeros(sample_count,model_feats*6);
    
    % Iterating over rows
    for i= 1:sample_count
       % Iterating over columns
        for j= 1:6
            sampled_ds(i,3*(j-1)+1) = mean(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
            sampled_ds(i,3*(j-1)+2) = rms(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)); 
            sampled_ds(i,3*(j-1)+3) = max(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)) - mean(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
        end
    end      

end
