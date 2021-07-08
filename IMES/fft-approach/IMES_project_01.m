%% Clear variables and windows
clear;
close all;
clc;

%% READ AND PLOT DATA

classf_feats = {'mcu_timestamp','bmi_accel_x', 'bmi_accel_y','bmi_accel_z','bmi_gyro_x','bmi_gyro_y','bmi_gyro_z'};
%ds_names = ["empty_ds", "wood_ds", "plastic_ds"];

empty_ds = read_and_clean("empty_imu.csv",classf_feats,30000,25000);
wood_ds = read_and_clean("wood_imu.csv",classf_feats,30000,25000);
plastic_ds = read_and_clean("plastic_imu.csv",classf_feats,30000,25000);

% Visualize data in time-domain
% Note - X-axis shows no. of points instead of time since all points are sampled at equal frequency i.e 1600 Hz
plot_ds(empty_ds);
plot_ds(wood_ds);
plot_ds(plastic_ds);

%% FOURIER TRANSFORM & DATA FILTERING
Fs = 1600;

%%%%%%%%%%%%%    
figure ;
L1 = size(empty_ds,1);
fft_empty_ax = fft(empty_ds(:,2));
P2_1 = abs(fft_empty_ax/L1);
P1_1 = P2_1(1:L1/2+1);
P1_1(2:end-1) = 2*P1_1(2:end-1);

f1 = Fs*(0:(L1/2))/L1;
plot(f1,P1_1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1_1(f)|')

%%%%%%%%%%%%%%
figure;
L2 = size(wood_ds,1);
fft_wood_ax = fft(wood_ds(:,2));
P2_2 = abs(fft_wood_ax/L2);
P1_2 = P2_2(1:L2/2+1);
P1_2(2:end-1) = 2*P1_2(2:end-1);

f2 = Fs*(0:(L2/2))/L2;
plot(f2,P1_2) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f2 (Hz)')
ylabel('|P1_2(f)|')

%%%%%%%%%%%%%%
figure;
L3 = size(plastic_ds,1);
fft_plastic_ax = fft(plastic_ds(:,2));
P2_3 = abs(fft_plastic_ax/L3);
P1_3 = P2_3(1:L3/2+1);
P1_3(2:end-1) = 2*P1_3(2:end-1);

f3 = Fs*(0:(L3/2))/L3;
plot(f3,P1_3) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f3 (Hz)')
ylabel('|P1_3(f)|')


%% PREPARE DATA FOR ML

% Convert dataset array to the corresponding feature matrix
empty_ds_ml = ds_with_features(empty_ds,3200,5);
wood_ds_ml = ds_with_features(wood_ds,3200,5);
plastic_ds_ml = ds_with_features(plastic_ds,3200,5);

% Stack all datasets together
all_ds_ml = [empty_ds_ml; wood_ds_ml; plastic_ds_ml];

% Initialize labels
labels = categorical([zeros(size(empty_ds_ml,1),1);ones(size(wood_ds_ml,1),1); 2*ones(size(plastic_ds_ml,1),1)]);

% Convet array to table 
ML_feat_table = array2table(all_ds_ml,'VariableNames',{'mean_ax','rms_ax','amp_ax','max_ax','min_ax','mean_ay','rms_ay','amp_ay','max_ay','min_ay','mean_az','rms_az','amp_az','max_az','min_az','mean_gx','rms_gx','amp_gx','max_gx','min_gx','mean_gy','rms_gy','amp_gy','max_gy','min_gy','mean_gz','rms_gz','amp_gz','max_gz','min_gz'});

% Add labels column
ML_feat_table.labels = labels;


%% APPENDIX: User-defined functions

% 1. read_and_clean : 
% WHAT IT DOES - returns the dataset as an array reading required features from csv file
% file_name : string, name of the csv file i.e. example.csv
% classf_feats : a vector containing required feature names as strings

function temp_ds_array = read_and_clean(file_name,classf_feats,start_pt,end_pt)
    
    % Import .csv files into datastore
    temp_ds = datastore(file_name,'TreatAsMissing','NA', 'MissingValue',0);
    
    % Convert selected features from datastore to array 
    temp_ds.SelectedVariableNames = classf_feats;
    temp_ds_array = table2array(readall(temp_ds));
    
    % Start timestamp counter from 0
    temp_ds_array(:,1) = temp_ds_array(:,1)- temp_ds_array(1,1);
    
    % Ignore rows with timestamp = 0
    % Why? - They do not have all the 1600 points for the 1st second (observed in all 3 datasets)
    %posn = temp_ds_array(:,1)==0;
    %temp_ds_array(posn,:)= [];
    temp_ds_array = temp_ds_array(start_pt:end-end_pt,:);
 
end

% 2. plot_ds :
% WHAT IT DOES - returns a 6 i.e. ax_ay,az,gx,gy,gz quantity sub-plot figure

function plot_ds(temp_ds)
    
    % Define labels for each sub-plot
    labels = ["accel_x", "accel_y","accel_z","gyro_x","gyro_y","gyro_z"];
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
            
            % Feature 1 - Mean
            sampled_ds(i,model_feats*(j-1)+1) = mean(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
            
            % Feature 2 - RMS
            sampled_ds(i,model_feats*(j-1)+2) = rms(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)); 
            
            % Feature 3 - Amplitude (Max value with respect to mean value)
            sampled_ds(i,model_feats*(j-1)+3) = max(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)) - mean(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));

            % Feature 4 - Max value
            sampled_ds(i,model_feats*(j-1)+4) = max(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
            
            % Feature 5 - Min value
            sampled_ds(i,model_feats*(j-1)+5) = min(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
            
            
        end
    end      
end
