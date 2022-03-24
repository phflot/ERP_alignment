function [tmp,tmp_original,v] = Iter_Al(inputImage,~,ref,iter_ref,sigma)

% if you want to plot time signals replace ~ by = time 

tic;            % start alignment
% Sigma2 = estimate_noise(mat2gray(inputImage));


%%  prefilter

% figure;imagesc(inputImage);colormap(hot);colorbar;   %check
% title('Input ABR matrix');
% ylabel('number of sweeps');
% xticklabels = time(1):2:time(end);
% xticks = linspace(1, size(time, 2), numel(xticklabels));
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% xlabel('time [ms]','FontSize', 15)
% yyaxis right;plot(mean(inputImage));

% BM3D
%==============================
%    [~, prefiltered] = BM3D(1, mat2gray(inputImage), Sigma, 0);
   [prefiltered] = BM3DSHARP(mat2gray(inputImage),sigma, 100, [], 0);


% figure;imagesc(prefiltered);colormap(gray);colorbar;  %check
% title('Prefilterd Matrix');
% ylabel('number of sweeps');
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% xlabel('time [ms]','FontSize', 15)

       
    %% 1 dimensional image alignment Philipp Flotho/David Thinnes 2019
  
    s_i = 10;   % 10
    it_i = 60;  % 60
    a_i = 1;  % 0.8
    as_i = 0.5;
    
    s_r = 10;    % 5
    it_r = 20; % 20
    a_r = 1;  % 0.8
    as_r = 0.5;  % 0.5
    
    
        % initial alignment to get the reference (select stable reference measurements):
    [~, v] = align_lines(...
        prefiltered, ...
        mat2gray(ref), ...      % prefiltered           %(40-60)
        'sigma',s_i, ...                                
        'iterations', it_i, ...                         % würde ich auf etwa 20-150 setzen
        'alpha', a_i, ...                               % Smoothness weight, falls der Code abbricht, einfach vergrößern, falls das Ergebnis zu gleichförmig ist einfach verkleinern
        'a_smooth', as_i);                              % 0.5 erlaubt Diskontinuitäten, 1 eher weiche Übergänge (brauchen aber meistens leider ein unteschiedliches alpha...)
    
    

    reg = horiz_alignment(inputImage, v);                      % inputImage
   
   [reg_filt, ~] = imgaussfiltaniso( mat2gray(reg), 4, 2, true);


% figure;imagesc(reg_filt);colormap(hot);colorbar;    %check
% title('Anistropic Gaussian');
% ylabel('number of sweeps');
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% xlabel('time [ms]','FontSize', 15)
% yyaxis right;plot(mean(reg_filt));



    
    for tt = 1:iter_ref
         
    % refinement:
    [~, v] = align_lines(...
        prefiltered, ...                                       % prefiltered    or registered ??     
        mat2gray(mean(reg_filt, 1)), ...
        'v', v, ...
        'sigma', s_r, ...
        'iterations', it_r, ...
        'alpha', a_r, ...
        'a_smooth', as_r);
    
  
    w(:, :, 1) = double(zeros(size(v)));
    w(:, :, 2) = v;
    
   
    
    registered = horiz_alignment(inputImage, v);                % inputImage    

    [reg_filt, ~] = imgaussfiltaniso( mat2gray(registered), 4, 2, true);
    end
    
     
   tmp_original = inputImage;
   tmp = registered;


 el_time = toc;
 fprintf(['elapsed_time Alignment= ',num2str(el_time),'\n']);
 
end

