%% Exercise 1: Ultrasound image display
% Demonstration of signal processing involved in ultrasound image display.
% This includes finding the enveloped through a Hilbert transform,
% compressing the data, and making the image interpolation. It also show
% how serveral frames can be combined into one movie.
% Exercise from http://bme.elektro.dtu.dk/31545/?exercises/exercise1/exercise1.html
% Assumes that the field_init routine has been run

%% Initilizing
close all, clear all, clc;
load rf_data_phantom;
RFdata = double(RFdata);

% Calculate the axial axis of the image [mm]
z = start_of_data + ((0:size(RFdata,1)-1).*c_sound)/(2*fs);

centerline=round(size(RFdata,2)/2);

%% Evelope

env = abs(hilbert(RFdata));

figure
hold on
plot(z*1000,RFdata(:,centerline))
plot(z*1000,env(:,centerline),'r')
title('Envelope using Hilberttransform')

%% Compression

env_dB = 20*log10(env/max(max(env))); % decibel scale
env_dB = env_dB-max(max(env_dB))+60; % shift
% env_dB = (env_dB+60).*(env_dB>-60)-60;
env_dB(env_dB<0) = 0; % cut off low values

figure
plot(z*1000,env_dB(:,centerline))
title('Log compressed, 60 dB dynamic range')

% env_gray = 127*(env_dB+60)/60;

%% Interpolation

% Setup parameters for make_tables
start_depth=start_of_data;
image_size=z(end)-start_of_data;
delta_r=c_sound/(2*fs);
N_samples=size(env,1);
theta_start=start_angle;
delta_theta=angle;
N_lines=size(env,2);
scaling=1;
Nz=N_samples;
Nx=N_lines;

% Calculate tables for interpolation
make_tables(start_depth, image_size, start_of_data, delta_r, ...
    N_samples, theta_start, delta_theta, N_lines, scaling, Nz, Nx);
% Calculate image
img_data = make_interpolation(uint32(env_dB));

%% Display

figure
imagesc([],z,img_data)
colormap(gray(128))
