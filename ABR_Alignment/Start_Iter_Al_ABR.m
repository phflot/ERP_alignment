
%% Start the variational_parataver for ABR

clc;
close all;
clear;
addpath(genpath('./'));

method = 'Mean';
% method = 'Woody';
% method = 'Cross';

% input
load('test60db_left.mat');
sweeps = inputImage(:,1:480);
load('time.mat');
Fs = 20000; 

% pad left and rigth array to avoid bounding effects
padsize = 20;
l = fliplr(sweeps(:,1:padsize));
r = fliplr(sweeps(:,size(sweeps,2)-(padsize-1):size(sweeps,2)));
sweeps = [l sweeps r];

time = time(1:480);  %24ms

  
 inputImage = sweeps;
 L1 = size(inputImage,1);
 L2 = size(inputImage,2);
 
 
 %% randomize the trial order
new_order = randperm(L1);
inputImage = inputImage(new_order,:);
 
%% define stable reference 
block = 50;
runs = floor(L1/block);
mean_inputs = [];
cc_scatter = [];


new_input = inputImage(1:block-1,:);    % set for first run
x1 =  0;
for ii = 1:runs-1
    switch method
        case {'Cross', 'cross'}
            ref  = CrossCorrRef(new_input,1);
            %         [ ref ] = RefCrossIter( inputImage );
            ref = new_input(ref,:);
            ref = mean(ref,1);
        case {'Mean' ,'mean'}
            ref = mean(new_input);
            %                 ref = T_testRef(inputImage);
        case {'Woody', 'woody'}
            ref = woody(new_input',[],[],'thornton','biased');
            ref = ref';
        otherwise
            ref = mean(new_input);
    end
    
    
    [tmp,tmp_original,v] = Iter_Al(new_input,time,ref,1,10);
    
    mean_inputs = [mean_inputs;mean(tmp)];
    
    
    if ii >1
        c = corrcoef(mean(inputImage),mean(mean_inputs));
        cc_value = c(1,2);
    else
        cc_value = 0;
    end
    cc_scatter = [cc_scatter; cc_value];
    
    x1 = x1+block;
    x2 = x1+(2*block);
    if ii < runs-1
        new_input = inputImage(x1:x2,:);
    end
    
    
    fprintf(['RUN:',num2str(ii),' Correlation with conv average = ',num2str(cc_value),'  '])
    
    
    if cc_value > 0.9
        break;
    end
    
end




[tmp_GA,tmp_original_GA,v_GA] = Iter_Al(mean_inputs,time,mean(mean_inputs),1,10);   

% random input trials for comparison
random_input_trials = inputImage(randperm(ii*block),:);

% delete the padding
inputImage_cropped = inputImage(:,padsize:end-(padsize+1));
mean_inputs_cropped = mean_inputs(:,padsize:end-(padsize+1));
tmp_GA_cropped = tmp_GA(:,padsize:end-(padsize+1));
v_GA_cropped = v_GA(:,padsize:end-(padsize+1));
random_input_trials_cropped = random_input_trials(:,padsize:end-(padsize+1));


figure;
subplot(1,3,1)
imagesc(inputImage_cropped);colormap(hot);colorbar;   %check
title('Input ABR matrix');
ylabel('number of sweeps');
xticklabels = time(1):2:time(end);
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 15)
yyaxis right;plot(mean(inputImage_cropped),'LineWidth',2);
ylim([min(mean(tmp_GA_cropped)) max(mean(tmp_GA_cropped))])

subplot(1,3,2)
imagesc(mean_inputs_cropped);colormap(hot);colorbar;   %check
title('Variatinal paratver ABRs');
ylabel('number of PAS');
xticklabels = time(1):2:time(end);
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 15)
yyaxis right;plot(mean(mean_inputs_cropped),'LineWidth',2);
ylim([min(mean(tmp_GA_cropped)) max(mean(tmp_GA_cropped))])

subplot(1,3,3)
imagesc(tmp_GA_cropped);colormap(hot);colorbar;   %check
title('Refined PAS');
ylabel('number of PAS');
xticklabels = time(1):2:time(end);
xticks = linspace(1, size(time, 2), numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
xlabel('time [ms]','FontSize', 15)
yyaxis right;plot(mean(tmp_GA_cropped),'LineWidth',2);
ylim([min(mean(tmp_GA_cropped)) max(mean(tmp_GA_cropped))])



figure;
subplot(1,2,2)
plot(time,mean(inputImage_cropped));
hold on;
plot(time,mean(mean_inputs_cropped));
plot(time,mean(tmp_GA_cropped,'omitnan'));
plot(time,mean(random_input_trials_cropped));
hold off;
title('Mean of all PAS blocks')
legend('mean input', 'mean PAS', 'mean refined PAS', 'random input trials')

subplot(1,2,1)
LL = ((1:ii)/runs)*100;
scatter(LL,cc_scatter,'+');
xlabel('% used input trials','FontSize', 15)
ylabel('pearsons correclation coeff.','FontSize', 15);
title('Correlation Aligned blocks and Avg input matrix','FontSize', 15)
lsline


c = corrcoef(mean(inputImage_cropped),mean(tmp_GA_cropped));
cc_end = c(1,2);
fprintf(['Correlation Aligned VS Input ABRs= ',num2str(cc_end)]);

