function [tmp,tmp_original,v] = Iter_ABR(inputImage,~,ref,iter_ref,sigma,denoising,half_filter)

% if you want to plot time signals replace ~ by = time 
tic;            % start alignment

%%  prefilter
% BM3D
%==============================
prefiltered = denoise('BM3D', inputImage, sigma, half_filter);
  
%% 1 dimensional image alignment Philipp Flotho/David Thinnes 2019
      
s_i = 15;
it_i = 80;
a_i = 0.8;
as_i = 0.5;

s_r = 10;
it_r = 50;
a_r = 0.4;
as_r = 0.5;

    
        % initial alignment to get the reference (select stable reference measurements):
    [~, v] = align_lines(...
        prefiltered, ...
        mat2gray(ref), ...      % prefiltered         
        'sigma',s_i, ...                                
        'iterations', it_i, ...                        
        'alpha', a_i, ...                             
        'a_smooth', as_i);                              
    

    registered = horiz_alignment(inputImage, v);                   
    initial_al = registered;


    %iterative application
    for tt = 1:iter_ref
    reg_filt = denoise(denoising, registered, sigma, half_filter);     
    % refinement:
    [~, v] = align_lines(...
        reg_filt, ...                                       % prefiltered    or registered ??     
        mat2gray(mean(reg_filt, 1)), ...
        'v', v, ...
        'sigma', s_r, ...
        'iterations', it_r, ...
        'alpha', a_r, ...
        'a_smooth', as_r);
    
  
    w(:, :, 1) = double(zeros(size(v)));
    w(:, :, 2) = v;
    
    registered = horiz_alignment(inputImage, v);                % inputImage    
    end
    
     
   tmp_original = inputImage;
   tmp = registered;


 el_time = toc;
 fprintf(['elapsed_time Alignment= ',num2str(el_time),'\n']);
 
end


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
