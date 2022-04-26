% Author: David Thinnes
% date:  04/22/2022
% Copyright 2022 by David Thinnes, All rights reserved.

% compare differnet prefilter algorithms with PSNR method
%% set path
clear
close all
clc

run('set_path_functionsERP');
load('ERPexample');

% generate synthetic ERP matrix
[synth_ERP, v_ref, reference, groundTruthMap] = generate_synth_ERP();

nr_iter = 1;
method = 'mean';
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


%% influence on iterative refinement
[Reg_aniso_1,V_aniso_1, initial_denoising_aniso,~] = Var_Alignment_comp_prefilter(input,method,1 ,Fs,'aniso',false);
[Reg_aniso_5,V_aniso_5, ~,~] = Var_Alignment_comp_prefilter(input,method,5 ,Fs,'aniso',false);

[Reg_BM3D_1,V_BM3D_1, initial_denoising_BM3D,initial_al] = Var_Alignment_comp_prefilter(input,method,1 ,Fs,'BM3D',false);
[Reg_BM3D_5,V_BM3D_5, ~,~] = Var_Alignment_comp_prefilter(input,method,5 ,Fs,'BM3D',false);

figure;
subplot(3,4,1)
imagesc(input);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('time [ms]','FontSize', 12);
ylabel('ERP trials','FontSize', 12)
title('noisy input')
subplot(3,4,2)
imagesc(initial_denoising_aniso)
title('proposed gaussian')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(3,4,3)
imagesc(Reg_aniso_1)
title('1 iteration');
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(3,4,4)
imagesc(Reg_aniso_5)
title('5 iterations')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(3,4,6)
imagesc(initial_denoising_BM3D)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
title('BM3D')
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(3,4,7)
imagesc(Reg_BM3D_1)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(3,4,8)
imagesc(Reg_BM3D_5)
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

subplot(3,4,11);
plot(time,mean(input));
hold on;
plot(time,mean(Reg_aniso_1)),
plot(time,mean(Reg_aniso_5));
xlim([time(1) time(end)])
grid on;
title('change in mean signal: BM3D')
legend('input','1 Iteration','5 Iterations','Location','NW');

subplot(3,4,12);
plot(time,mean(input));
hold on;
plot(time,mean(Reg_BM3D_1)),
plot(time,mean(Reg_BM3D_5));
xlim([time(1) time(end)])
grid on;
title('change in mean signal: gaussian')
legend('input','1 Iteration','5 Iterations','Location','NW');


