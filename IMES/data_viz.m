%%
clear;
close all;
clc;

%%



%%
function importdata
empty_imu_ds = datastore('empty_imu.csv','TreatAsMissing','NA','MissingValue',0);

classification_features = {'mcu_timestamp', 'bmi_gyro_y'};
empty_imu_ds.SelectedVariableNames = classification_features;
empty_ds_clean = readall(empty_imu_ds)

%%
function plot_features
plot(empty_ds_clean.mcu_timestamp, empty_ds_clean.bmi_gyro_y )

%% 


%% 