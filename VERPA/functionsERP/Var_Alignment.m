function [Reg,V,initial_denoising, initial_al] = Var_Alignment(inputImage,method,iter_ref,Fs,denoising,half_filter,plot_results)
% Author: David Thinnes
% date:  03/27/2023
% Copyright 2023 by David Thinnes, All rights reserved.

% description:
% the function corrects the latency shift of evoked potentials by choosing a stable
% reference and realigning the 2D Sweep Image of the data

% it performs 1 dimensional alignment according to the selected reference
% method: cross correlation, mean of all image rows 
% default: 'Mean'

% input: 
%===============================================================================================
% type :        input Image (Sweep matrix, Single trial matrix)
% method :      'Mean', 'Cross', 'woody', 'thornton'   
% denoising:    'ansio','aspr','diffusion','BM3D','wavelet','nlm', 'regu'
% iter_ref :    iterative refinement
% Fs :          Sampling Frequency of the signal

% output
%===============================================================================================
% Reg =         the realigned image after image registration
%   V =         the displacement field of all pixels, samples

% please refere to:
%===============================================================================================
% variational_aligner Philipp Flotho, 2020
% Philipp Flotho David Thinnes, F. I. Corona Strauss, J. F. Vibell  & D.J. Strauss    Fast Variational Method for the Estimation and Quantization of non–flat Displacements in 1D Signals with Applications in Neuroimaging, 2020
% David Thinnes Philipp Flotho, F. I. Corona Strauss, D. J. Strauss  & J.F. Vibell   Compensation of P300 Latency Jitter using fast variational 1D Displacement Estimation, 2020


%% define image size

 % size of inputImage
 L1 = size(inputImage,1);
 L2 = size(inputImage,2);

%% create time vector
L = (L2-1)/Fs;
time  = 0:1/Fs:L;
time = time.*1000;      % in ms
time = time -100;       % first 100 ms are baseline
    
tic;            % start alignment program

%% define the reference
switch method
    
    case 'Cross', case 'cross'
        [ ref ] = CrossCorrRef(inputImage,1);
        ref = inputImage(ref,:);
        ref = mean(ref,1);
    case 'Mean' , case 'mean'
        ref = mean(inputImage);
    case {'Woody','woody'}
        ref = woody(inputImage',[],[],'woody','unbiased');
        ref = ref';
    case {'Thornton','thornton'}
        ref = woody(inputImage',[],[],'thornton','biased');
        ref = ref';

    otherwise
        ref = mean(inputImage); % default
end



%% define parameters for the alignment program
s_i = 15;
it_i = 80;
a_i = 0.8;
as_i = 0.5;

s_r = 10;
it_r = 50;
a_r = 0.4;
as_r = 0.5;

reg_filt = denoise(denoising, inputImage, 5, half_filter);


initial_denoising = reg_filt;


% initial alignment to get the reference (select stable reference measurements):
[~, v] = align_lines(...
    reg_filt, ...
    mat2gray(ref), ...
    'sigma',s_i, ...
    'iterations', it_i, ...                         % about 20-150
    'alpha', a_i, ...                               
    'a_smooth', as_i);                              
registered = horiz_alignment(inputImage, v);               % inputImage

initial_al = registered;



% iterative application possible
for tt = 1:iter_ref
    reg_filt = denoise(denoising, registered, 2, half_filter);
   
    % refinement:
    [~, v_tmp] = align_lines(...
        reg_filt, ...
        mat2gray(mean(reg_filt, 1)), ...
        ...'v', v, ...
        'sigma', s_r, ...
        'iterations', it_r, ...
        'alpha', a_r, ...
        'a_smooth', as_r);
    
    
    v = v + v_tmp;

    w(:, :, 1) = double(zeros(size(v)));
    w(:, :, 2) = v;
    

    registered = horiz_alignment(inputImage, v);          % inputImage
    

end

%% Save the results
tmp_original = inputImage;
tmp = registered;

Reg = registered;
V = v;

% show elapsed time of the program
elapsedtime = toc;
fprintf(['Elapsed Time Alignment =' num2str(elapsedtime) 'sec \n']);

%% Plot the results
if plot_results == true
    
    figure('units','normalized','outerposition',[0 0 1 1]);         % fullsize figure
ROI1 = 250;
ROI2 = 400;
x = [ROI1 ROI2  ROI2  ROI1];
y = [(min(mean(tmp))-0.2) (min(mean(tmp))-0.2) (max(mean(tmp))+0.2) (max(mean(tmp))+0.2)];

subplot(2,3,1)
imagesc(tmp_original);
xticklabels = time(1):100:600;
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('time [ms]','FontSize',15);
ylabel('ERP trials','FontSize', 15)
title('Pre Alignment','FontSize', 15);
colorbar;

subplot(2,3,2)
imagesc(initial_denoising);
xticklabels = time(1):100:600;
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('time [ms]','FontSize',15);
ylabel('ERP trials','FontSize', 15)
title(['Denoised ERPs: ',denoising],'FontSize', 15);
colorbar;


subplot(2,3,3);
imagesc(tmp);
title('Post Alignment','FontSize', 15);
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('time [ms]','FontSize',15);
ylabel('ERP trials','FontSize', 15)
colorbar;

subplot(2,3,4);
patch(x,y, [0.85 0.85 0.85]);
hold on;
title('Average of all trials','FontSize', 15);
plot(time,mean(inputImage),'LineWidth', 2); 
mean_tmp = mean(tmp, 1, 'omitnan');
plot(time,mean_tmp,'LineWidth', 2);
legend('ROI','Pre-Alignment','Post-Alignment','Location','NW');
xlabel('time [ms]','FontSize', 15)
grid on;
hold off;

% plot displacements for a certain time frame (over a certain trace)
subplot(2,3,[5 6]);

Dis = squeeze(mean(v(:, find(single(time) == ROI1):find(single(time) == ROI2) ), 2));
plot(Dis, 'LineWidth', 2);
xlim([1 L1]);
xlabel('number of Sweeps','FontSize', 15)
ylabel('displacements','FontSize', 15)
title(['Mean displacements between  ', num2str(ROI1), ':', num2str(ROI2),'ms'],'FontSize', 15)
grid on;

end

clear w v registered
end


% select the denoising filter
function  reg_filt = denoise(denoising, inputImage, sigma, half_filter)
if nargin < 4
    half_filter = false;
end
switch denoising
    case 'aniso'
        % filter with anisotropic gaussian
        if half_filter == true
        [reg_filt, ~] = imgaussfiltaniso( inputImage, sigma, 4 * sigma, half_filter);
        else
        [reg_filt, ~] = imgaussfiltaniso( inputImage, sigma, 4 * sigma);
        end
    case 'aspr'
        % filter with amplitude and phase regularization
        reg_filt = ASPRdenoising(inputImage);
    case 'diffusion'
        % filter with anisotropic diffusion
        options = get_aniso_coh_diff_options();
        options.diffusivity = 'weickert';
        reg_filt = aniso_coh_diff_filt(inputImage, options);
        figure;imagesc(reg_filt);

    case 'wavelet'
        % filter with wavelet denoising
        reg_filt = wdenoise2(inputImage);
    case 'BM3D'
        % filter with block matching 3D denoising
%         sigma = estimate_noise(inputImage);
         sigma = 10;
         [~, reg_filt] = BM3D(1, mat2gray(inputImage), sigma);


    case 'nlm'
        % filter with non local means
        [reg_filt,~] = imnlmfilt(inputImage);
    case 'regu'
        % filter with anisotropic phase-map denoising
        reg_filt = regu(inputImage,2,25,0.01);
end

end
