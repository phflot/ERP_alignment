% Author: David Thinnes
% date:  04/22/2022
% Copyright 2022 by David Thinnes, All rights reserved.

% compare differnet prefilter algorithms with PSNR method
%% set path
clear
close all
clc

run('set_path_functions');
load('ERPexample');

% generate synthetic ERP matrix
[synth_ERP, v_ref, reference, groundTruthMap] = generate_synth_ERP();

nr_iter_low = 1;
nr_iter_high = 4;
method = 'mean';
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


%% influence on iterative refinement
[Reg_aniso_1,V_aniso_1, initial_denoising_aniso,~] = Var_Alignment(input,method,nr_iter_low ,Fs,'aniso',false,false);
[Reg_aniso_5,V_aniso_5, ~,~] = Var_Alignment(input,method,nr_iter_high,Fs,'aniso',false,false);

[Reg_BM3D_1,V_BM3D_1, initial_denoising_BM3D,initial_al] = Var_Alignment(input,method,nr_iter_low ,Fs,'BM3D',true,false);
[Reg_BM3D_5,V_BM3D_5, ~,~] = Var_Alignment(input,method,nr_iter_high,Fs,'BM3D',true,false);

figure;
subplot(3,4,1)
imagesc(input);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('time [ms]','FontSize', 15);
ylabel('ERP trials','FontSize', 15)
title('noisy input','FontSize', 15)
xtickangle(45)
subplot(3,4,2)
imagesc(initial_denoising_aniso)
title('proposed gaussian','FontSize', 15)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(3,4,3)
imagesc(Reg_aniso_1)
title('1 iteration','FontSize', 15);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(3,4,4)
imagesc(Reg_aniso_5)
title('5 iterations','FontSize', 15)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(3,4,6)
imagesc(initial_denoising_BM3D)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
title('BM3D','FontSize', 15)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(3,4,7)
imagesc(Reg_BM3D_1)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(3,4,8)
imagesc(Reg_BM3D_5)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)

subplot(3,4,11);
plot(time,mean(input),'LineWidth',1);
hold on;
plot(time,mean(Reg_aniso_1),'LineWidth',1),
plot(time,mean(Reg_aniso_5),'LineWidth',1);
xlim([time(1) time(end)])
grid on;
title('BM3D','FontSize', 15)
xlabel('time [ms]','FontSize',15)
legend('input','1 Iteration',[num2str(nr_iter_high),'iteration'],'Location','NW','FontSize', 10);
ylim([-1 1.5]);

subplot(3,4,12);
plot(time,mean(input),'LineWidth',1);
hold on;
plot(time,mean(Reg_BM3D_1),'LineWidth',1),
plot(time,mean(Reg_BM3D_5),'LineWidth',1);
xlim([time(1) time(end)])
grid on;
title('Aniso gaussian','FontSize', 15)
xlabel('time [ms]','FontSize',15)
% legend('input','1 Iteration',[num2str(nr_iter_high),'iteration'],'Location','NW','FontSize', 10);
ylim([-1 1.5]);

