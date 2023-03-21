% Author: David Thinnes
% date:  04/22/2022
% Copyright 2022 by David Thinnes, All rights reserved.

%% set path
clear
close all
clc

run('set_path_functions');
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

ROI1 = 250;
ROI2 = 400;
%% influence on iterative refinement
[Reg_BM3D_1,V_BM3D_1, initial_denoising_BM3D,initial_al] = Var_Alignment(input,method,1 ,Fs,'aniso',false,false);

Dis = squeeze(mean(V_BM3D_1(:, find(single(time) == ROI1):find(single(time) == ROI2) ), 2));
Dis_ms = Dis*5;

figure;
subplot(4,1,1)
imagesc(groundTruthMap);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('time [ms]','FontSize', 12);
ylabel('shifted ERP trials','FontSize', 12)
title('GroundTruth map')
h = colorbar;
title(h, 'norm. potential')


subplot(4,1,2)
imagesc(synth_ERP);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
ylabel('shifted ERP trials','FontSize', 12)
title('Pre-Alignment (noise added)')
h = colorbar;
title(h, 'potential')
hold on;
yyaxis right;
plot(mean(synth_ERP),'LineWidth',2,'Color',[128,128,128]/255)
ax = gca;
ax.YAxis(2).Color = [128,128,128]/255;
ylim([-1 +1.5]);


subplot(4,1,3)
imagesc(V_BM3D_1);
ylabel('1D displacements','FontSize', 12)
title('Calculated displacements')
h = colorbar;
title(h, 'Samples')
x = 90 +Dis;
x_loc =x/Fs;
y = 1:1:length(Dis);
hold on;
plot(x,y,'k','LineWidth',1);




subplot(4,1,4)
imagesc(Reg_BM3D_1);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
ylabel('aligned ERP trials','FontSize', 12)
title('Post-Alignment')
h = colorbar;
title(h, 'potential')
hold on;
yyaxis right;
plot(mean(synth_ERP),'-','LineWidth',2,'Color',[128,128,128]/255)
plot(mean(Reg_BM3D_1),'-','LineWidth',2,'Color','r')
ax = gca;
ax.YAxis(2).Color = [128,128,128]/255;
ylim([-1 +1.5]);
