%% Start the Alignment on the whole ABR Image

clc;
close all;
clear;

addpath(genpath('./'));

method = 'Mean';
% method = 'Woody';


% input
load('0db.mat');
sweeps = sweeps.Left(:,1:480);
load('time.mat');
Fs = 20000; 


padsize = 20;
l = fliplr(sweeps(:,1:padsize));
r = fliplr(sweeps(:,size(sweeps,2)-(padsize-1):size(sweeps,2)));
sweeps = [l sweeps r];

time = time(1:480);  %24ms


m =  sweeps;
snum = 10;


% partaver
[row,column]=size(m);
pval=round(row/snum)-1;
sweeps=zeros(pval,column);
for k=1:1:pval
   sweeps(k,:)=mean(m(((k-1)*snum+1):(k*snum),:));
end

% % artefact
% TFACT= 15; %for artefact removal
% [sweeps] = artefact(sweeps,TFACT);
% sweeps = - sweeps; % change polarity (Cz = ref)



Mean_sweeps = mean(sweeps);
  

 inputImage = sweeps;
 L1 = size(inputImage,1);
 L2 = size(inputImage,2);

% %% create time vector
% L = (L2-1)/Fs;
% time  = 0:1/Fs:L;
% time = time.*1000;      % in ms#


% prefilter anisotropic gaussian
% Sigma2 = estimate_noise(mat2gray(inputImage));

tic;            % start alignment
 
%% define stable reference 

switch method
    
    case {'Cross','cross'}
        [ ref ] = CrossCorrRef(inputImage,1);
%         [ ref ] = RefCrossIter( inputImage );
        ref = inputImage(ref,:);
        ref = mean(ref,1);

    case {'Mean', 'mean'}
        ref = mean(inputImage);
%                 ref = T_testRef(inputImage);
    case {'Woody','woody'}
        ref = woody(inputImage',[],[],'thornton','biased');
        ref = ref';

    otherwise
        ref = inputImage(1,:);
end




%%  prefilter

figure;imagesc(inputImage);colormap(hot);colorbar;   %check
title('Input ABR matrix');
ylabel('number of sweeps');
xticklabels = time(1):2:time(end);
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 15)
yyaxis right;plot(mean(inputImage));

% BM3D
%==============================
%    [~, prefiltered] = BM3D(1, mat2gray(inputImage), Sigma, 0);

[prefiltered] = BM3DSHARP(mat2gray(inputImage),0.1, 500, [], 0);
   

% 
% figure;imagesc(prefiltered);colormap(gray);colorbar;  %check
% title('Prefilterd Matrix');
% ylabel('number of sweeps');
% set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
% xlabel('time [ms]','FontSize', 15)


    s_i = 10;   % 10
    it_i = 80;  % 60
    a_i = 1;  % 0.8
    as_i = 0.5;
    
    s_r = 5;    % 5
    it_r = 20; % 20
    a_r = 1;  % 0.8
    as_r = 0.5;  % 0.5

    
    
    %% 1 dimensional image registration program Philipp Flotho 2019
    
        % initial alignment to get the reference (select stable reference measurements):
    [~, v] = align_lines(...
        prefiltered, ...
        mat2gray(ref), ...      % prefiltered           %(40-60)
        'sigma',s_i, ...                                
        'iterations', it_i, ...                         % würde ich auf etwa 20-150 setzen
        'alpha', a_i, ...                               % Smoothness weight, falls der Code abbricht, einfach vergrößern, falls das Ergebnis zu gleichförmig ist einfach verkleinern
        'a_smooth', as_i);                              % 0.5 erlaubt Diskontinuitäten, 1 eher weiche Übergänge (brauchen aber meistens leider ein unteschiedliches alpha...)
    
    reg = horiz_alignment(inputImage, v);                      % inputImage
    

    [reg_filt, ~] = imgaussfiltaniso( reg, 4, 5, true);


    iter_ref = 1;

    
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

    [reg_filt, ~] = imgaussfiltaniso( registered, 4, 5, true);
    end
    
     
   tmp_original = inputImage;
   tmp = registered;
  

   
% delete the padding
registered_cropped = registered(:,padsize:end-(padsize+1));
tmp_original_cropped = inputImage(:,padsize:end-(padsize+1));
reg_cropped = reg(:,padsize:end-(padsize+1));


v_cropped = v(:,padsize:end-(padsize+1));
   
   
   
figure;imagesc(registered_cropped);colormap(hot);colorbar;    %check
title('Aligned ABR post partaver');
ylabel('number of sweeps');
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 15)

yyaxis right;plot(mean(registered));

      %% plot the results

figure('units','normalized','outerposition',[0 0 1 1]);         % fullsize figure
colormap(hot);
subplot(2,3,1)
imagesc(tmp_original_cropped);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
  ax = gca;
  ax.FontSize = 12;
xlabel('time [ms]','FontSize', 21)
ylabel('number of Sweeps','FontSize', 21)
title('Input ABR','FontSize', 22);
%     colorbar;
    
subplot(2,3,2);
imagesc(reg_cropped);
  ax = gca;
  ax.FontSize = 12;
title('Initial Alignment','FontSize', 22);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 21)
% ylabel('number of Sweeps','FontSize', 15)

subplot(2,3,3);
imagesc(tmp);
  ax = gca;
  ax.FontSize = 12;
title('Refinement','FontSize', 22);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 22)
% ylabel('number of Sweeps','FontSize', 15)
%     colorbar;
%     subplot(2,2,[3,4]);
subplot(2,3,[4 5 6]);

plot(time,Mean_sweeps(:,padsize:end-(padsize+1)),'LineWidth',3,'Color',[0.5 0.5 0.5]);
hold on;
plot(time,mean(reg(:,padsize:end-(padsize+1))),'LineWidth',2.5,'Color','b','LineStyle',':');
hold on;                                 % Set Zeros To ‘NaN’
mean_tmp = mean(tmp(:,padsize:end-(padsize+1)));
plot(time,mean_tmp,'LineWidth',2.5,'Color','r','LineStyle','--')
      ax = gca;
  ax.FontSize = 12;
title('Mean Pre-Alignment vs Mean Post-Alignment','FontSize', 22)
    
    xlim([time(1) time(end)]);
    xlabel('time [ms]','FontSize', 21)
    ylabel('potential [\muV]','FontSize', 21)
%     legend('mean input','mean output');
    legend({'mean input','mean intital','mean refinement'},'FontSize', 12);
    grid on;
    hold off;
   
   
figure;
subplot(1,2,1);
imagesc(v_cropped);
title('Displacement map')
ylabel('sweeps')
xlabel('Samples')
subplot(1,2,2);
dis_med = median(v_cropped(:,80:200))/(Fs/1000);
dis_mean  = mean(v_cropped(:,80:200))/(Fs/1000);
time_int = time(80:200);
plot(time_int,dis_med);
hold on;
plot(time_int,dis_mean);
plot(time_int(22), dis_med(22), 'o')
plot(time_int(22), dis_mean(22), 'o')
xlabel('time interval [ms]')
ylabel('time alignment of sweeps [ms]')
title('Alignment of matrix')
grid on;
legend('median','mean')
    

    