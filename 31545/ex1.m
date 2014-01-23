%% Exercise 1: Ultrasound image display
% Demonstration of signal processing involved in ultrasound image display. This includes finding the enveloped through a Hilbert transform, compressing the data, and making the image interpolation. It also show how serveral frames can be combined into one movie.
% Exercise from http://bme.elektro.dtu.dk/31545/?exercises/exercise1/exercise1.html
% Assumes that the field_init routine has been run

clc; close all; clear all;

%% Initilizing
load rf_data_phantom;
RFdata = double(RFdata);

depth = start_of_data + ((0:size(RFdata,1)-1)*c_sound)/(2*fs);

%% Evelope and compression

env = abs(hilbert(RFdata));

figure
hold on
plot(depth*1000,RFdata(:,140))
plot(depth*1000,env(:,140),'r')
title('Envelope using Hilberttransform')

env_dB = 20*log10(env/max(max(env)));
env_dB = (env_dB+60).*(env_dB>-60)-60;

figure
plot(depth*1000,env_dB(:,140))

env_grey = 127*(env_dB+60)/60;

%% Interpolation

image_size = [depth(end),depth(end)];
delta_r = c_sound/(2*fs);
make_tables(start_of_data,image_size,...
    start_of_data,delta_r,size(RFdata,1),...
    start_angle,angle,size(RFdata,2),...
    1,400,400);
img_data = make_interpolation(uint32(env_dB));

imagesc(img_data,gray(128))
