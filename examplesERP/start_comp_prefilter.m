% Author: David Thinnes
% date:  03/03/2022
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

Fs = 200;
input = synth_ERP;
sigma = 5;
half_filter = false;


% compare the different prefilter algorithms 
initial_denoising_aniso = denoise('aniso', input, sigma, half_filter);
initial_denoising_aspr = denoise('aspr', input, sigma, half_filter);
initial_denoising_diffusion = denoise('diffusion', input, sigma, half_filter);
initial_denoising_wavelet = denoise('wavelet', input, sigma, half_filter);
initial_denoising_BM3D = denoise('BM3D', input, sigma, half_filter);
initial_denoising_nlm = denoise('nlm', input, sigma, half_filter);
initial_denoising_regu = denoise('regu', input, sigma, half_filter);


%% create time vector
L1 = size(input,1);
L2 = size(input,2);
L = (L2-1)/Fs;
time  = 0:1/Fs:L;
time = time.*1000;      % in ms
time = time -100;       % first 100 ms are baseline

xticklabels = time(1):100:600;
xticks = linspace(1, size(time, 2), numel(xticklabels));

figure;
subplot(2,4,1)
imagesc(input)
PSNRdb = PSNR(groundTruthMap,input); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['input: PSNR ',num2str(PSNRdb), 'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,2);
imagesc(initial_denoising_aniso)
PSNRdb = PSNR(groundTruthMap,initial_denoising_aniso); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['aniso gauss: PSNR ',num2str(PSNRdb), 'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,4);
imagesc(initial_denoising_aspr)
PSNRdb = PSNR(groundTruthMap,initial_denoising_aspr); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['aspr: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,5);
imagesc(initial_denoising_diffusion),
PSNRdb = PSNR(groundTruthMap,initial_denoising_diffusion); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['diffusion: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,8);
imagesc(initial_denoising_wavelet)
PSNRdb = PSNR(groundTruthMap,initial_denoising_wavelet); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['wavelet: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,7);
imagesc(initial_denoising_nlm)
PSNRdb = PSNR(groundTruthMap,initial_denoising_nlm); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['nlm: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,3);
imagesc(initial_denoising_regu)
PSNRdb = PSNR(groundTruthMap,initial_denoising_regu); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['regu: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xtickangle(45)
subplot(2,4,6);
imagesc(initial_denoising_BM3D)
PSNRdb = PSNR(groundTruthMap,initial_denoising_BM3D); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));
title(['BM3D: PSNR ',num2str(PSNRdb),'dB'])
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
sgtitle('pre-alignment denosing')
xtickangle(45)




function  filt = denoise(denoising, inputImage, sigma, half_filter)
if nargin < 4
    half_filter = false;
end
switch denoising
    case 'aniso'
        % filter with anisotropic gaussian
        if half_filter == true
            [filt, ~] = imgaussfiltaniso( inputImage, sigma, 4 * sigma, half_filter);
        else
            [filt, ~] = imgaussfiltaniso( inputImage, sigma, 4 * sigma);
        end
    case 'aspr'
        % filter with amplitude and phase regularization
        filt = ASPRdenoising(inputImage);
    case 'diffusion'
        % filter with anisotropic diffusion
        options = get_aniso_coh_diff_options();
        options.diffusivity = 'weickert';
        filt = aniso_coh_diff_filt(inputImage, options);

    case 'wavelet'
        % filter with wavelet denoising
        filt = wdenoise2(inputImage);
    case 'BM3D'
        % filter with block matching 3D denoising
        %         sigma = estimate_noise(inputImage);
        sigma = 10;
        [~, filt] = BM3D(1, mat2gray(inputImage), sigma);

    case 'nlm'
        % filter with non local means
        [filt,~] = imnlmfilt(inputImage);
    case 'regu'
        % filter with anisotropic phase-map denoising
        filt = regu(inputImage,2,25,0.01);
end
end


