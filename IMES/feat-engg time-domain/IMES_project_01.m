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
%plot_ds(empty_ds);
%plot_ds(wood_ds);
%plot_ds(plastic_ds);

%% PREPARE DATA FOR ML

model_feats = 3; 

% Convert dataset array to the corresponding feature matrix
empty_ds_ml = ds_with_features(empty_ds,3200,model_feats);
wood_ds_ml = ds_with_features(wood_ds,3200,model_feats);
plastic_ds_ml = ds_with_features(plastic_ds,3200,model_feats);

% Stack all datasets together
all_ds_ml = [empty_ds_ml; wood_ds_ml; plastic_ds_ml];

% Initialize labels
labels = categorical([zeros(size(empty_ds_ml,1),1);ones(size(wood_ds_ml,1),1); 2*ones(size(plastic_ds_ml,1),1)]);
featNames = ["mean_ax","mean_ay","mean_az","mean_gx","mean_gy","mean_gz","rms_ax","rms_ay","rms_az","rms_gx","rms_gy","rms_gz","amp_ax","amp_ay","amp_az","amp_gx","amp_gy","amp_gz","max_ax","max_ay","max_az","max_gx","max_gy","max_gz","min_ax","min_ay","min_az","min_gx","min_gy","min_gz","mode_ax","mode_ay","mode_az","mode_gx","mode_gy","mode_gz"];

% Index iterator over the featNames cell
indx = model_feats*6; 
avail_feats = featNames(1:indx);

% Convet array to table 
ML_feat_table = array2table(all_ds_ml,'VariableNames',avail_feats);

% Add labels column
ML_feat_table.labels = labels;

%% TRAIN & TEST MODEL

% Split into train test dataset
c = cvpartition(ML_feat_table.labels,'HoldOut',0.25); %10% use for testing data
triidx = training(c);
testingdata = ML_feat_table(~triidx,:);
training_validation_data = ML_feat_table(triidx,:);

%Use classification learning apps
%Select training_validation_data to train and validate

%% Export Model

%Classify the testing dataset using trainedModel
pred = trainedModel.predictFcn(testingdata);

accuracy = (1-mean(pred~=testingdata.labels))*100;
fprintf("Test accuracy is %f ", accuracy)
%% APPENDIX: User-defined functions

% 1. read_and_clean : 
% WHAT IT DOES - returns the dataset as an array reading required features from csv file
% file_name : string, name of the csv file i.e. example.csv
% classf_feats : a vector containing required feature names as strings
% start_pts : number of points to drop from the start of the data
% end_pts : number of points to drop from the end of the data

function temp_ds_array = read_and_clean(file_name,classf_feats,start_pts,end_pts)
    
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
    
    % Drop points during transient state of the machine i.e. during the
    % start and the stopppag
    temp_ds_array = temp_ds_array(start_pts:end-end_pts,:);
 
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
            
            if (model_feats<7 && model_feats>0)
                
                count = 0;
                
                if (count<model_feats)
                    % Feature 1 - Mean
                    sampled_ds(i,model_feats*(j-1)+1) = mean(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)); 
                    count = count + 1;
                end
                
                if (count < model_feats)
                    % Feature 2 - RMS
                    sampled_ds(i,model_feats*(j-1)+2) = rms(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)); 
                    count = count +1;
                end
                
                if (count < model_feats)
                    % Feature 3 - Amplitude (Max value with respect to mean value)
                    sampled_ds(i,model_feats*(j-1)+3) = max(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j)) - mean(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
                    count = count +1;
                end
                
                if (count < model_feats)
                    % Feature 4 - Max
                    sampled_ds(i,model_feats*(j-1)+4) = max(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
                    count = count +1;
                end

                if (count < model_feats)
                    % Feature 5 - Min
                    sampled_ds(i,model_feats*(j-1)+5) = min(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
                    count = count +1;
                end
                
                if (count < model_feats)
                    % Feature 6 - Mode
                    sampled_ds(i,model_feats*(j-1)+6) = min(temp_ds((i-1)*sampling_factor+1:i*sampling_factor,j));
                    count = count +1;
                end
            
                %disp(count)    
                
            else
                disp("INVALID entry: Number of features either exceeds or falls short of no. of min. features ")
            end
            
        end
    end      
end
