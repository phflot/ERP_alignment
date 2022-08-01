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

% compare the different prefilter algorithms 
[Reg_aniso,V_aniso, initial_denoising_aniso,initial_al] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'aniso',true);
[Reg_aspr,V_aspr, initial_denoising_aspr,~] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'aspr',true);
[Reg_diffusion,V_diffusion,initial_denoising_diffusion,~] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'diffusion',true);
[Reg_donoho,V_donoho, initial_denoising_donoho,~] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'wavelet',true);
[Reg_BM3D,V_BM3D, initial_denoising_BM3D,~] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'BM3D',true);
[Reg_nlm,V_nlm, initial_denoising_nlm,~] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'nlm',true);
[Reg_regu,V_regu, initial_denoising_regu,~] = Var_Alignment_comp_prefilter(input,method,nr_iter ,Fs,'regu',true);





%% create time vector
L1 = size(input,1);
L2 = size(input,2);
L = (L2-1)/Fs;
time  = 0:1/Fs:L;
time = time.*1000;      % in ms
time = time -100;       % first 100 ms are baseline

xticklabels = time(1):100:time(end);
xticks = linspace(1, size(time, 2), numel(xticklabels));

figure;
subplot(2,4,1)
imagesc(input)
PSNRdb = PSNR(groundTruthMap,input); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['noisy input: PSNR ',num2str(PSNRdb), 'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,2);
imagesc(initial_denoising_aniso)
PSNRdb = PSNR(groundTruthMap,initial_denoising_aniso); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['aniso gauss: PSNR ',num2str(PSNRdb), 'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,4);
imagesc(initial_denoising_aspr)
PSNRdb = PSNR(groundTruthMap,initial_denoising_aspr); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['aspr: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,5);
imagesc(initial_denoising_diffusion),
PSNRdb = PSNR(groundTruthMap,initial_denoising_diffusion); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['diffusion: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,8);
imagesc(initial_denoising_donoho)
PSNRdb = PSNR(groundTruthMap,initial_denoising_donoho); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['wavelet: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,7);
imagesc(initial_denoising_nlm)
PSNRdb = PSNR(groundTruthMap,initial_denoising_nlm); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['nlm: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,3);
imagesc(initial_denoising_regu)
PSNRdb = PSNR(groundTruthMap,initial_denoising_regu); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['regu: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
subplot(2,4,6);
imagesc(initial_denoising_BM3D)
PSNRdb = PSNR(groundTruthMap,initial_denoising_BM3D); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['BM3D: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
sgtitle('pre-alignment denosing')



