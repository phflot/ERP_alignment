% MATLAB code for implementation of iterative adaptive image deblurring
% algorithm:
% A. Kheradmand and P. Milanfar, “A general framework for regularized, similarity-based image restoration,” 
% IEEE Transactions on Image Processing, vol. 23, no. 12, pp. 5136–5151, Dec 2014.
% This is experimental software. It is provided for non-commercial research purposes only. Use at your own risk. No warranty is implied by 
% this distribution. Copyright © 2014 by University of California.

%% producing the blurred image
%%cunstructing the blur kernel (PSF)
clc; clear all;
addpath('BM3D');

scenario =1;
switch scenario
    case 1  %Gaussian blur
        width = 25;
        h1=fspecial('gaussian',width,1.6); 
        sigma = sqrt(1); blurcase = 'sym';
    case 2 %motion blur
        load mblur; h1 = double(mblur);
        h1 = h1/sum(h1(:)); 
         sigma = sqrt(1);
        blurcase = 'motion';
    case 3 %out-of-focus blur
        width = 7; 
        h1 = fspecial('disk', width); sigma = sqrt(1);
        blurcase = 'sym';
end

z = double(imread('motobikes.png'));
N = size(z,1); M = size(z,2); Hf = psf2otf(h1, [N M]);

y_blurred = zeros(N, M);
for i=1:3
    y_blurred(:,:,i) = real(ifft2(fft2(z(:,:,i)).*Hf));
end

randn('state',1);

y = y_blurred + sigma * randn(size(y_blurred)); %Blurred noisy image

% Noise variance estimation
sigma_hatt = zeros(1,3);
for i = 1:3
    sigma_hatt(i) = estimate_noise(y(:,:,i));
end

%% kernel coefficient computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_deblur = zeros(size(y));
nOut = 3; %Number of outer itertions
for kk = 1:3
    fprintf('Initial denoising color channel %0.0f ... \n ',kk);
    [PSNR1, y_est] = BM3D(y_blurred(:,:,kk)/255, y(:,:,kk)/255, sigma_hatt(kk), 0);
    zh = y_est * 255; close all; z1 = zh; % zh is initial denoised image
    fprintf('Iterations color channel %0.0f ... \n ',kk);
    for jj = 1:nOut
        if (jj~=1)
            z1 = z_final;
        end
        
        switch blurcase
            case 'sym'
                if sigma == sqrt(0.2)
          numiter = 100 - 30*(jj-1); eta = 0.003; beta = 0.2; hr = 5.5; %parameters for sigma^2 = 0.2
                else
                numiter = 100 - 30*(jj-1); eta = 0.008; beta = 0.001; hr= 7.5; %parameters for sigma^2 = 1
                end
            case 'motion'
                if sigma == sqrt(0.2)
               numiter = 100 - 30*(jj-1); eta = 0.006; beta = 0.4; hr = 6; %parameters for sigma^2 = 0.2
                else numiter = 80 - 30*(jj-1); eta = 0.01; beta = 0.01; hr = 6.5; %parameters for sigma^2 = 1
                end
        end
        wrad = 5; ksz= 5; % Parameters for NLM kernel similarity weights
        [K_rec] = Integral_NLM(z1,(ksz-1)/2,wrad,hr); % Computing the kernel similarity matrix

        %% Constructing filterin matrix W %%%%%%%%%%
        y1 = y(:,:,kk);
        [x,res] = scaling_sp(K_rec,1e-1); % Fast Sinkhorn Algorithm to get diag(C^{-1/2})
        Kernell = sparse(1:N*M, 1:N*M, x, N*M, N*M, numel(x)); clear x; % Diagonal matrix C^{-1/2}
        K_rec = Kernell*K_rec*Kernell; % Creating Filtering matrix 
        clear Kernell;
        y_estt = reshape(K_rec*y1(:),N,M);
        z_est = zh; % Initilaiztion of the inner CG iterations

        if (jj~=1)
            z_est = z_final;
        end

%% CG iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [PMSE, z_final] = kernel_CG_W (y1,z_est,...
            Hf, K_rec, beta, eta, numiter, zh);
        clear z_est;
    end
    z_deblur(:,:,kk) = z_final;
    fprintf('---------------------------------------------- \n');
end

% Displaying the results
figure; imshow(uint8(z),[0 255]); title('Original Image');
figure; imshow(uint8(y),[0 255]); title('blurred noisy input');
figure; imshow(uint8(z_deblur), [0 255]); title('deblurred image');
% Quantitative evaluation of the algorithm
PSNRdb = PSNR(z, z_deblur); fprintf('PSNR of the output image is %f\n',mean(PSNRdb));