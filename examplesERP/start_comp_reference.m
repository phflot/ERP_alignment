% Author: David Thinnes
% date:  08/02/2022
% Copyright 2022 by David Thinnes, All rights reserved.

% compare differnet references (mean, woody thornton)
%% set path
clear
close all
clc

run('set_path_functionsERP');
load('ERPexample');

% generate synthetic ERP matrix
[synth_ERP, v_ref, reference, groundTruthMap] = generate_synth_ERP();

nr_iter = 1;
Fs = 200;
input = synth_ERP;


%% compare the denoising methods
% create time vector
L1 = size(input,1);
L2 = size(input,2);
L = (L2-1)/Fs;
time  = 0:1/Fs:L;
time = time.*1000;      % in ms
time = time -100;       % first 100 ms are baseline

xticklabels = time(1):100:time(end);
xticks = linspace(1, size(time, 2), numel(xticklabels));

% references
ref_mean = mean(input);

ref_woody = woody(input',[],[],'woody','unbiased');
ref_woody = ref_woody';

ref_thornton = woody(input',[],[],'thornton','biased');
ref_thornton = ref_thornton';


%% influence of the choice of reference
method = 'mean';
[Reg_mean,V_mean, initial_denoising_BM3D,initial_al] = Var_Alignment_comp_prefilter(input,method,1 ,Fs,'BM3D',false,false);
method = 'woody';
[Reg_woody,V_woody, ~,~] = Var_Alignment_comp_prefilter(input,method,1 ,Fs,'BM3D',false,false);
method = 'thornton';
[Reg_thornton,V_thornton, ~,~] = Var_Alignment_comp_prefilter(input,method,1 ,Fs,'BM3D',false,false);

figure;
subplot(3, 2, 1);
plot(ref_mean,'color','g','LineWidth',2);
title('Template for initial alignment','FontSize',15);
ylabel('norm. amplitude','FontSize',15);
grid on;
legend('mean');


subplot(3, 2, 3);
plot(ref_woody,'color','r','LineWidth',2);
ylabel('norm. amplitude','FontSize',15);
grid on;
legend('woody');

subplot(3, 2, 5);
plot(ref_thornton,'color','k','LineWidth',2);
ylabel('norm. amplitude','FontSize',15);
grid on;
legend('thornton');
xlabel('time [ms]','FontSize',15);

subplot(3, 2, [2 4 6])
plot(mean(input),'color','b','LineWidth',2);
hold on;
plot(mean(Reg_mean),'color','g','LineWidth',2);
plot(mean(Reg_woody),'color','r','LineWidth',2);
plot(mean(Reg_thornton),'color','k','LineWidth',2);
title('Post Alignment','FontSize',15);
ylabel('norm. amplitude','FontSize',15);
xlabel('time [ms]','FontSize',15);
grid on;
legend({'conventional','mean','woody','thornton'});