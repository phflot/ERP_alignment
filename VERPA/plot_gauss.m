% Author: David Thinnes
% date:  05/02/2022
% Copyright 2022 by David Thinnes, All rights reserved.

% plot the anisotropic gaussian

%% set path
clear
close all
clc

run('set_path_functionsERP');
load('ERPexample');

input = ERPexample;

[~, G] = imgaussfiltaniso( input, 5, 15, false);

figure;
subplot(121)
imagesc(G);
colorbar;
xlabel('size','FontSize',12)
ylabel('size','FontSize',12)
title('full','FontSize',12)

[~, G] = imgaussfiltaniso( input, 5, 15, true);

subplot(122)
imagesc(G)
colorbar;
xlabel('size','FontSize', 12)
ylabel('size','FontSize',12)
title('half','FontSize',12)

sgtitle('Angular weighted gaussian','FontSize',12);
