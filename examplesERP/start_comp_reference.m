% Author: David Thinnes
% date:  03/03/2023
% Copyright 2022 by David Thinnes, All rights reserved.

%% compare differnet references (mean, woody thornton)
%% set path
clear
% close all
clc

run('set_path_functions');
load('ERPexample');

% generate synthetic ERP matrix
[synth_ERP, v_ref, reference, groundTruthMap] = generate_synth_ERP();

nr_iter = 1;
Fs = 200;
input = synth_ERP;


% create time vector
L1 = size(input,1);
L2 = size(input,2);
L = (L2-1)/Fs;
time  = 0:1/Fs:L;
time = time.*1000;      % in ms
time = time -100;       % first 100 ms are baseline

xticklabels = time(1):100:600;
xticks = linspace(1, size(time, 2), numel(xticklabels));

% references
ref_mean = mean(input);

ref_woody = woody(input',[],[],'woody','unbiased');
ref_woody = ref_woody';

ref_thornton = woody(input',[],[],'thornton','biased');
ref_thornton = ref_thornton';


%% influence of the choice of reference
method = 'mean';
[Reg_mean,V_mean, initial_denoising_BM3D,initial_al] = Var_Alignment(input,method, nr_iter ,Fs,'aniso',false,false);
method = 'woody';
[Reg_woody,V_woody, ~,~] = Var_Alignment(input,method,nr_iter ,Fs,'aniso',false,false);
method = 'thornton';
[Reg_thornton,V_thornton, ~,~] = Var_Alignment(input,method,nr_iter,Fs,'aniso',false,false);

figure;
subplot(3, 2, 1);
plot(time,ref_mean,'color','g','LineWidth',2);
title('Template for initial alignment','FontSize',15);
% ylabel('norm. amplitude','FontSize',15);
grid on;
ylim([-1 1.5])
legend('mean','location','nw');


subplot(3, 2, 3);
plot(time,ref_woody,'color','k','LineWidth',2);
% ylabel('norm. amplitude','FontSize',15);
grid on;
ylim([-1 1.5])
legend('woody','location','nw');

subplot(3, 2, 5);
plot(time,ref_thornton,'color','r','LineWidth',2);
% ylabel('norm. amplitude','FontSize',15);
grid on;
ylim([-1 1.5])
legend('thornton','location','nw');
xlabel('time [ms]','FontSize',15);

subplot(3, 2, [2 4 6])
plot(time,mean(input),'color','b','LineWidth',2);
hold on;
plot(time,mean(Reg_mean),'color','g','LineWidth',2);
plot(time,mean(Reg_woody),'color','k','LineWidth',2);
plot(time,mean(Reg_thornton),'color','r','LineWidth',2);
title('Post Alignment','FontSize',15);
ylabel('norm. amplitude','FontSize',15);
xlabel('time [ms]','FontSize',15);
grid on;
legend({'conventional','mean','woody','thornton'},'location','nw');


figure;
subplot(1,2,1)
imagesc(synth_ERP)
title('Synthetic ERP with jitter','FontSize',15);
ylabel('ERP trials','FontSize',15);
xlabel('time [ms]','FontSize',15);
% xticklabels = time(1):100:600;
% xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

subplot(1,2,2)
plot(time,reference,':','color','k','LineWidth',2);
hold on;
plot(time,mean(input),'color','b','LineWidth',2);
plot(time,mean(Reg_mean),'color','g','LineWidth',2);
plot(time,ref_woody,'color','k','LineWidth',2);
plot(time,ref_thornton,'color','r','LineWidth',2);
title('Cross-Correlation VS VM Alignment','FontSize',15);
ylabel('norm. amplitude','FontSize',15);
xlabel('time [ms]','FontSize',15);
grid on;
legend({'GroundTruth','conventional','our method','woody','thornton'},'location','nw');


%% Align SNR
figure;

reg_mean_align = align_CC(mean(Reg_mean, 1), reference);
reg_thorn_align = align_CC(ref_thornton, reference);
reg_woody_align = align_CC(ref_woody, reference);
input_mean_align = align_CC(mean(input, 1), reference);

ref = reference(:, 15:end-15);

disp(psnr(input_mean_align, ref))
disp(psnr(reg_woody_align, ref))
disp(psnr(reg_thorn_align, ref))
disp(psnr(reg_mean_align, ref))


hold on
plot(ref,':','color','k','LineWidth',1);
plot(input_mean_align,'color','b','LineWidth',1);
plot(reg_mean_align,'color','g','LineWidth',1);
plot(reg_woody_align,'color','k','LineWidth',1);
plot(reg_thorn_align,'color','r','LineWidth',1);

title('SNR Alignment','FontSize',15);
grid on;
legend({'GroundTruth','conventional','our method','woody','thornton'},'location','nw');

function f = align_CC(f, ref)
    [R, lag] = xcorr(f, ref);
    [~, idx] = max(R);
    f = circshift(f, -lag(idx));
    f = f(:, 15:end-15);
end