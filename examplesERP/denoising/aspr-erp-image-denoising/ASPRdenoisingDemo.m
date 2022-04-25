
% Filename : ASPRdenoisingDemo.m
% Date     : 15.02.2018
% Authors  : Manuel C. Kohl

clc;
clear;
close all;
addpath(genpath('./'));

% Loading ERP image

load('ERPimage.mat');
[nResponses, nSamples] = size(ERPimage);

% Denoising ERP image

denoisedERPimage = ASPRdenoising(ERPimage);

% Plotting data

figure;
figurePosition = get(gcf, 'position');
figurePosition(3) = 2 * figurePosition(3);
set(gcf, 'position', figurePosition);
centerfig;
subplot(1, 2, 1);
imagescnan(responseTime, 1:nResponses, ERPimage);
colormap(linearInferno());
colorBar = colorbar();
zLabel = get(colorBar, 'xlabel');
set(zLabel, 'String', 'Amplitude U(n, t) [\muV]');
xlabel('Time t [ms]');
ylabel('Responses n');
title('Unprocessed ERP image');
subplot(1, 2, 2);
imagescnan(responseTime, 1:nResponses, denoisedERPimage);
colormap(linearInferno());
colorBar = colorbar();
zLabel = get(colorBar, 'xlabel');
set(zLabel, 'String', 'Amplitude U(n, t) [\muV]');
xlabel('Time t [ms]');
ylabel('Responses n');
title('Denoised ERP image');
